/*
 * lovo.go, part of gochem.
 *
 *
 * Copyright 2021 Raul Mera rauldotmeraatusachdotcl
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 *
 *
 */

package align

import (
	//	"bufio"

	"fmt"
	"log"
	"math"
	"runtime"
	"sort"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/traj/dcd"
	"github.com/rmera/gochem/traj/xtc"
	v3 "github.com/rmera/gochem/v3"
)

//ConcTraje is an interface for a trajectory that can be read concurrently.
type ConcTrajCloser interface {
	chem.ConcTraj
	//closes the trajectory
	Close()
}

//Options contains various options for the LOVOnMostRigid and RMSDTraj functions
type Options struct {
	Begin int
	Skip  int
	Cpus  int
	//The following are ignored by the RMSDTraj function
	AtomNames    []string
	Chains       []string
	NMostRigid   int
	WriteTraj    string  //the name of a file to where the aligned trajectory will be written. Nothing will be written if empty.
	LessThanRMSD float64 //instead of using a N for the most rigid, selects all residues with RMSD (the square root of the RMSD) < LessThanRMSD, in A. If >0, this overrrides NMostRigid
	MinimumN     int     //The smallest N we are willing to have. Only valid if LessThanRMSD is in use.
}

//DefaultOptions return reasonable options for atomistic trajectories.
//It prepares a superposition of alpha carbons (CA) with all logical CPUs,
//trying to use for the superposition all CAs with RMSD lower than 1.0 A.
func DefaultOptions() *Options {
	r := new(Options)
	r.Begin = 0
	r.Skip = 0
	r.Cpus = runtime.NumCPU()
	//all available CPUs
	r.AtomNames = []string{"CA"}
	r.Chains = nil
	r.NMostRigid = -1
	r.LessThanRMSD = 1.0
	r.MinimumN = 10 //just a reasonable value.
	return r
}

//DefaultCGOptions returns reasonable options for Martini trajectories.
func DefaultCGOptions() *Options {
	r := DefaultOptions()
	r.AtomNames = []string{"BB"}
	r.Skip = 1000 //seems reasonable for CG, those trajectories are _long_.
	return r

}

//Sets O.N to represent the perc percent of the residues
//in the sequence. It also requires seqlen, the total number of
//residues in the system
func (O *Options) SetRigidPercent(perc int, seqlen int) {
	frac := float64(perc) / 100
	O.NMostRigid = int(frac * float64(seqlen))

}

//LOVOReturn contains the information returned by LOVOnMostRigid
type LOVOReturn struct {
	N          int
	FramesRead int       //how many frames where actually read from the trajectory
	Natoms     []int     //IDs of the N most rigid atoms
	Nmols      []int     //IDs for the N most rigid residues/molecules
	RMSD       []float64 //the RMSD for all residues, not only the N most rigid.
	Iterations int       //The iterations that were needed for convergency.

}

//String returns a string representation of the LOVOReturn object.
func (L *LOVOReturn) String() string {
	return fmt.Sprintf("N: %d, Natoms: %v, Nmols: %v, RMSD: %v, Frames read: %d, Iterations needed: %d", L.N, L.Natoms, L.Nmols, L.RMSD, L.FramesRead, L.Iterations)
}

//PyMOLSel returns a string of text to create a PyMOL
//selection with the L.N most rigid residues from a LOVO
//calculation.
func (L *LOVOReturn) PyMOLSel() string {
	pymolsele := "select rigid,"
	for i, v := range L.Nmols {
		pymolsele += fmt.Sprintf(" resi %d ", v)
		if i < len(L.Nmols)-1 {
			pymolsele += "or"
		}
	}
	return pymolsele

}

func opentraj(name string) (ConcTrajCloser, error) {
	t := strings.Split(name, ".")
	format := t[len(t)-1]
	var err error
	var ret ConcTrajCloser
	switch strings.ToLower(format) {
	case "xtc":
		ret, err = xtc.New(name)
	case "dcd":
		ret, err = dcd.New(name)
	default:
		ret = nil
		err = fmt.Errorf("Trajectory fromat not supported. Must have concurrent read support in goChem")
	}
	return ret, err

}

//tolerance is the minimum amount of residues we are willing to use in the superposition.
func mobileLimit(limit float64, tolerance int, RMSD *rMSD) (int, float64) {
	var n = 0
	for n < tolerance {
		n = 0
		for ; n < RMSD.Len(); n++ {
			if RMSD.RMSD(n) >= limit {
				n--
				break
			}
		}
		limit *= 1.1
	}
	return n, limit
}

//LOVO returns slices with the indexes of the n atoms and residues
//most rigid from mol, along the traj trajectory, considering only atoms with names in atomnames (normally, PDB names are used)
//and belonging to one of the chains in chains (only the first slice given in chains will be considered, if nothing is given,
//atoms from any chain will be considred.
//If you use this function in your research, please cite the reference for the LOVO alignment method:
//10.1371/journal.pone.0119264.
//If you use this function in a to-consumer program, please print a message asking to cite the reference when the function is used.
func LOVO(mol chem.Atomer, ref *v3.Matrix, trajname string, o *Options) (*LOVOReturn, error) {
	//we first do one iteration aligning the whole thing.
	//these are atom, not residue indexes!
	const defaultTopRigid int = 10 //by default we obtain the 1/topdefaultTopRigid top rigid
	var fullindexes []int = res2atoms(mol, o.AtomNames, o.Chains, nil)
	var printtraj string
	printtraj = o.WriteTraj
	o.WriteTraj = ""
	indexesold := make([]int, len(fullindexes))
	copy(indexesold, fullindexes)
	var RMSD *rMSD
	var n int = len(fullindexes)
	var totalframes int
	var itercount int
	var disagreement int
	var prevlim float64
	var rmsd []float64
	for {
		traj, err := opentraj(trajname)
		if err != nil {
			return nil, err
		}
		if disagreement == 1 && printtraj != "" {
			o.WriteTraj = printtraj
		}
		rmsd, totalframes, err = RMSDTraj(mol, ref, traj, indexesold[0:n], o)
		if err != nil {
			return nil, err
		}
		traj.Close()
		itercount++
		RMSD = newRMSD(fullindexes, rmsdFilter(rmsd, fullindexes))
		RMSD.SortBy("rmsd")
		var newlim float64
		if o.LessThanRMSD > 0 {
			n, newlim = mobileLimit(o.LessThanRMSD, o.MinimumN, RMSD)
		}
		indexes := RMSD.MolIDsCopy()
		if sameElementsInt(indexes[0:n], indexesold[0:n]) && (newlim == o.LessThanRMSD || approxEq(newlim, prevlim, 0.01)) {
			break //converged
		}
		prevlim = newlim
		disagreement = disagreementInt(indexes[0:n], indexesold[0:n])
		indexesold = indexes
	}
	//Now indexes old should countain what we want
	RMSD.SortBy("molid")
	r := &LOVOReturn{N: n, Natoms: indexesold[0:n], Nmols: iDs2MolIDs(mol, indexesold[0:n]), RMSD: RMSD.rMSDs, FramesRead: totalframes, Iterations: itercount}
	return r, nil

}

type rmsdandcoords struct {
	rmsd   []float64
	coords *v3.Matrix
}

func concproc(v chan *v3.Matrix, ref, t *v3.Matrix, r chan *rmsdandcoords, indexes []int) {
	//the errors that Super and MemRMSD can return all
	//would be programming erros in this case, so I'll just turn them into panics.
	var err2 error
	c := <-v
	cr, err2 := chem.Super(c, ref, indexes, indexes)
	if err2 != nil {
		panic(err2.Error())
	}
	rmsd, err2 := chem.MemPerAtomRMSD(cr, ref, nil, nil, t)
	if err2 != nil {
		panic(err2.Error() + fmt.Sprintf("Lengths: %d", t.NVecs()))
	}
	ret := &rmsdandcoords{rmsd: rmsd, coords: cr}
	r <- ret

}

//RMSDTraj returns the RMSD for all atoms in a structure, averaged over the trajectory traj.
//Frames are read in sets of cpus at the time: The read starts at the set begin/cpus and we skip a set every skip/cpu
//(the numbers are, of course, rounded). The RMSDs and the total number of frames read are returned, together with an error or nil.
func RMSDTraj(mol chem.Atomer, ref *v3.Matrix, traj chem.ConcTraj, indexes []int, o *Options) ([]float64, int, error) {
	RMSD := make([]float64, mol.Len())
	var err error
	chunksize := o.Cpus
	chunk := make([]*v3.Matrix, chunksize)
	tmps := make([]*v3.Matrix, chunksize)
	for i, _ := range chunk {
		chunk[i] = v3.Zeros(mol.Len())
		tmps[i] = v3.Zeros(mol.Len())
	}
	var wtraj *dcd.DCDWObj
	if o.WriteTraj != "" {
		wtraj, err = dcd.NewWriter(o.WriteTraj, mol.Len())
		if err != nil {
			log.Printf("Couldn't open %s for writing. Will not write aligned trajectory", o.WriteTraj) //no error handling here. If we can't open the file, we just don't write anything.
		} else {
			defer wtraj.Close()

		}

	}
	begin := o.Begin / chunksize //how many chunks to skip before begining.
	skip := o.Skip / chunksize   //how many chunks to skip
	results := make([]chan *rmsdandcoords, chunksize)
	for i, _ := range results {
		results[i] = make(chan *rmsdandcoords)
	}
	var chans []chan *v3.Matrix
	var chunksread int = 0
	for read := 0; ; read++ {
		if (read) < begin || (read)%(skip+1) != 0 {
			c := make([]*v3.Matrix, chunksize) //all nil
			_, err = traj.NextConc(c)
			if err != nil {
				break
			}
		}
		chans, err = traj.NextConc(chunk)
		if err != nil {
			break
		}
		chunksread++
		for i, v := range chans {
			go concproc(v, ref, tmps[i], results[i], indexes)

		}
		for _, res := range results {
			trmsd := <-res
			for i, _ := range RMSD {
				RMSD[i] += trmsd.rmsd[i]
			}
			if wtraj != nil {
				wtraj.WNext(trmsd.coords)
			}
		}

	}
	if err == nil {
		log.Println("Finished reading trajectory without a EOF signal. This shouldn't happen.")
	} else if _, ok := err.(chem.LastFrameError); !ok {
		return nil, -1, fmt.Errorf("Error while reading trajectory, read %d sets of frames. Error: %s", chunksread, err.Error())
	}
	totalframes := float64(chunksread * chunksize)
	for i, v := range RMSD {
		RMSD[i] = v / totalframes
	}
	return RMSD, int(totalframes), nil
}

type rMSD struct {
	//Note that the default methods and basis vary with each program, and even
	//for a given program they are NOT considered part of the API, so they can always change.
	//This is unavoidable, as methods change with time
	molIDs       []int
	rMSDs        []float64
	residues     int
	lessbyMolIDs func(i, j int) bool
	lessbyRMSDs  func(i, j int) bool
	sorting      string
}

func newRMSD(residues []int, iniRMSD ...[]float64) *rMSD {
	ret := new(rMSD)
	ret.residues = len(residues)
	ret.molIDs = make([]int, len(residues))
	copy(ret.molIDs, residues)
	if len(iniRMSD) > 0 {
		M := iniRMSD[0]
		if len(M) != ret.residues {
			panic("Inconsistent data, if given, the number of initial RMSD values must match the number of residues")
		}
		ret.rMSDs = make([]float64, len(M))
		copy(ret.rMSDs, M)

	}
	ret.lessbyMolIDs = func(i, j int) bool { return ret.molIDs[i] < ret.molIDs[j] }
	ret.lessbyRMSDs = func(i, j int) bool { return ret.rMSDs[i] < ret.rMSDs[j] }
	return ret
}

//MolID returns the molecule identifier of the nth residue of the structure.
func (m *rMSD) MolIDsCopy(rets ...[]int) []int {
	var ret []int
	if len(rets) != 0 {
		ret = rets[0]
	} else {

		ret = make([]int, len(m.molIDs))
	}
	copy(ret, m.molIDs)
	//	//By hand
	//	for i, v := range m.molIDs {
	//		ret[i] = v
	//
	//	}
	return ret //I'll just let it panic if the number given is too large.
}

func (m *rMSD) Swap(i, j int) {
	m.molIDs[i], m.molIDs[j] = m.molIDs[j], m.molIDs[i]
	m.rMSDs[i], m.rMSDs[j] = m.rMSDs[j], m.rMSDs[i]
}

//returns the ith RMSD in the current order
func (m *rMSD) RMSD(i int) float64 {
	return m.rMSDs[i]
}

func (m *rMSD) Len() int {
	return m.residues
}

func (m *rMSD) Less(i, j int) bool {
	switch m.sorting {
	case "rmsd":
		return m.lessbyRMSDs(i, j)
	default:
		return m.lessbyMolIDs(i, j)

	}

}

func (m *rMSD) SortBy(sorting string) {
	m.sorting = strings.ToLower(sorting)
	sort.Stable(m)
}

//helper functions

func approxEq(i, j, epsilon float64) bool {
	if math.Abs(i)-math.Abs(j) <= epsilon {
		return true
	}
	return false

}

//returns true if t1 and t2 have the same elements
//(whether or not in the same order) and false otherwise.
func sameElementsInt(t1, t2 []int) bool {
	if len(t1) != len(t2) {
		return false
	}
	for _, v := range t1 {
		if !isInInt(v, t2) {
			return false
		}
	}
	return true

}

func disagreementInt(t1, t2 []int) int {
	if len(t1) != len(t2) {
		log.Printf("inCommonInt: Slices differ in size! %d %d", len(t1), len(t2))
		return -1
	}
	dis := make([]int, 0, len(t1))
	for _, v := range t1 {
		if !isInInt(v, t2) {
			dis = append(dis, v)
		}
	}
	return len(dis)
}

//isInInt is a helper for the RamaList function,
//returns true if test is in container, false otherwise.
func isInInt(test int, container []int) bool {
	if container == nil {
		return false
	}
	for _, i := range container {
		if test == i {
			return true
		}
	}
	return false
}

//Same as the previous, but with strings.
func isInString(test string, container []string) bool {
	if container == nil {
		return false
	}
	for _, i := range container {
		if test == i {
			return true
		}
	}
	return false
}

//Some internal convenience functions.

//retuns the indexes in mol of the atoms with name in names and residue id (MolID) in res.
//if res is empty returns all atoms with name in names.
func res2atoms(mol chem.Atomer, names, chains []string, res []int) []int {
	ret := make([]int, 0, len(res))
	for i := 0; i < mol.Len(); i++ {
		a := mol.Atom(i)
		if (res == nil || isInInt(a.MolID, res)) && isInString(a.Name, names) {
			if len(chains) == 0 || isInString(a.Chain, chains) {
				ret = append(ret, i)
			}
		}
	}
	return ret

}

//returns a slice of RMSDs where the values with an index not in indexes have been removed.
func rmsdFilter(rmsd []float64, indexes []int) []float64 {
	ret := make([]float64, 0, len(indexes))
	for i, v := range rmsd {
		if isInInt(i, indexes) {
			ret = append(ret, v)
		}
	}
	return ret

}

func iDs2MolIDs(mol chem.Atomer, indexes []int) []int {
	ret := make([]int, 0, len(indexes))
	for _, v := range indexes {
		at := mol.Atom(v)
		if !isInInt(at.MolID, ret) {
			ret = append(ret, at.MolID)
		}
	}
	return ret

}
