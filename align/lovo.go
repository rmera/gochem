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
	"runtime"
	"sort"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/dcd"
	v3 "github.com/rmera/gochem/v3"
	"github.com/rmera/gochem/xtc"
)

//ConcTraje is an interface for a trajectory that can be read concurrently.
type ConcTrajCloser interface {
	chem.ConcTraj
	//closes the trajectory
	Close()
}

//Options contains various options for the LOVOnMostRigid and MSDTraj functions
type Options struct {
	Begin int
	Skip  int
	Cpus  int
	//The following are ignored by the MSDTraj function
	AtomNames  []string
	Chains     []string
	NMostRigid int
	WriteTraj  string //the name of a file to where the aligned trajectory will be written. Nothing will be written if empty.
}

//DefaultOptions return reasonable options for atomistic trajectories.
func DefaultOptions() *Options {
	r := new(Options)
	r.Begin = 0
	r.Skip = 0
	r.Cpus = runtime.NumCPU()
	//all available CPUs
	r.AtomNames = []string{"CA"}
	r.Chains = nil
	r.NMostRigid = -1
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
	MSD        []float64 //the MSD for all residues, not only the N most rigid.
	Iterations int       //The iterations that were needed for convergency.

}

//String returns a string representation of the LOVOReturn object.
func (L *LOVOReturn) String() string {
	return fmt.Sprintf("N: %d, Natoms: %v, Nmols: %v, MSD: %v, Frames read: %d, Iterations needed: %d", L.N, L.Natoms, L.Nmols, L.MSD, L.FramesRead, L.Iterations)
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
	}
	return ret, err

}

//LOVOnMostRigid returns slices with the indexes of the n atoms and residues
//most rigid from mol, along the traj trajectory, considering only atoms with names in atomnames (normally, PDB names are used)
//and belonging to one of the chains in chains (only the first slice given in chains will be considered, if nothing is given,
//atoms from any chain will be considred.
//If you use this function in your research, please cite the reference for the LOVO alignment method:
//10.1371/journal.pone.0119264.
//If you use this function in a to-consumer program, please print a message asking to cite the reference when the function is used.
func LOVOnMostRigid(mol chem.Atomer, ref *v3.Matrix, trajname string, o *Options) (*LOVOReturn, error) {
	//we first do one iteration aligning the whole thing.
	//these are atom, not residue indexes!
	const defaultTopRigid int = 10 //by default we obtain the 1/topdefaultTopRigid top rigid
	var fullindexes []int = res2atoms(mol, nil, o.AtomNames, o.Chains)
	var printtraj string
	printtraj = o.WriteTraj
	o.WriteTraj = ""
	traj, err := opentraj(trajname)
	if err != nil {
		return nil, err
	}
	msd, _, err := MSDTraj(mol, ref, traj, fullindexes, o)
	if err != nil {
		return nil, err
	}
	MSD := newMSD(fullindexes, msdFilter(msd, fullindexes))
	MSD.SortBy("msd")

	indexesold := MSD.MolIDsCopy()
	n := o.NMostRigid
	seqlen := mol.Atom(mol.Len() - 1).MolID //this is just a crude approximation for lack of something better. It will fail if the sequence is divided in different chains.
	if n <= 0 {
		n = seqlen / defaultTopRigid
	}
	var totalframes int
	var itercount int
	var disagreement int
	for {
		traj, err := opentraj(trajname)
		if err != nil {
			return nil, err
		}
		if disagreement == 1 && printtraj != "" {
			o.WriteTraj = printtraj
		}
		msd, totalframes, err = MSDTraj(mol, ref, traj, indexesold[0:n], o)
		if err != nil {
			return nil, err
		}
		itercount++
		MSD := newMSD(fullindexes, msdFilter(msd, fullindexes))
		MSD.SortBy("msd")
		indexes := MSD.MolIDsCopy()
		if sameElementsInt(indexes[0:n], indexesold[0:n]) {
			break //converged
		}
		disagreement = disagreementInt(indexes[0:n], indexesold[0:n])
		indexesold = indexes
	}
	//Now indexes old should countain what we want
	MSD.SortBy("molid")
	r := &LOVOReturn{N: n, Natoms: indexesold[0:n], Nmols: iDs2MolIDs(mol, indexesold[0:n]), MSD: MSD.mSDs, FramesRead: totalframes, Iterations: itercount}
	return r, nil

}

type msdandcoords struct {
	msd    []float64
	coords *v3.Matrix
}

func concproc(v chan *v3.Matrix, ref, t *v3.Matrix, r chan *msdandcoords, indexes []int) {
	//the errors that Super and MemMSD can return all
	//would be programming erros in this case, so I'll just turn them into panics.
	var err2 error
	c := <-v
	cr, err2 := chem.Super(c, ref, indexes, indexes)
	if err2 != nil {
		panic(err2.Error())
	}
	msd, err2 := chem.MemMSD(cr, ref, nil, nil, t)
	if err2 != nil {
		panic(err2.Error() + fmt.Sprintf("Lengths: %d", t.NVecs()))
	}
	ret := &msdandcoords{msd: msd, coords: cr}
	r <- ret

}

//MSDTraj returns the MSD for all atoms in a structure, averaged over the trajectory traj.
//Frames are read in sets of cpus at the time: The read starts at the set begin/cpus and we skip a set every skip/cpu
//(the numbers are, of course, rounded). The MSDs and the total number of frames read are returned, together with an error or nil.
func MSDTraj(mol chem.Atomer, ref *v3.Matrix, traj chem.ConcTraj, indexes []int, o *Options) ([]float64, int, error) {
	MSD := make([]float64, mol.Len())
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
	results := make([]chan *msdandcoords, chunksize)
	for i, _ := range results {
		results[i] = make(chan *msdandcoords)
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
			tmsd := <-res
			for i, _ := range MSD {
				MSD[i] += tmsd.msd[i]
			}
			if wtraj != nil {
				wtraj.WNext(tmsd.coords)
			}
		}

	}
	if err == nil {
		log.Println("Finished reading trajectory without a EOF signal. This shouldn't happen.")
	} else if _, ok := err.(chem.LastFrameError); !ok {
		return nil, -1, fmt.Errorf("Error while reading trajectory, read %d sets of frames. Error: %s", chunksread, err.Error())
	}
	totalframes := float64(chunksread * chunksize)
	for i, v := range MSD {
		MSD[i] += v / totalframes
	}
	return MSD, int(totalframes), nil
}

type mSD struct {
	//Note that the default methods and basis vary with each program, and even
	//for a given program they are NOT considered part of the API, so they can always change.
	//This is unavoidable, as methods change with time
	molIDs       []int
	mSDs         []float64
	residues     int
	lessbyMolIDs func(i, j int) bool
	lessbyMSDs   func(i, j int) bool
	sorting      string
}

func newMSD(residues []int, iniMSD ...[]float64) *mSD {
	ret := new(mSD)
	ret.residues = len(residues)
	ret.molIDs = make([]int, len(residues))
	copy(ret.molIDs, residues)
	if len(iniMSD) > 0 {
		M := iniMSD[0]
		if len(M) != ret.residues {
			panic("Inconsistent data, if given, the number of initial MSD values must match the number of residues")
		}
		ret.mSDs = M

	}
	ret.lessbyMolIDs = func(i, j int) bool { return ret.molIDs[i] < ret.molIDs[j] }
	ret.lessbyMSDs = func(i, j int) bool { return ret.mSDs[i] < ret.mSDs[j] }
	return ret
}

//MolID returns the molecule identifier of the nth residue of the structure.
func (m *mSD) MolIDsCopy(rets ...[]int) []int {
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

func (m *mSD) Swap(i, j int) {
	m.molIDs[i], m.molIDs[j] = m.molIDs[j], m.molIDs[i]
	m.mSDs[i], m.mSDs[j] = m.mSDs[j], m.mSDs[i]
}

func (m *mSD) Len() int {
	return m.residues
}

func (m *mSD) Less(i, j int) bool {
	switch m.sorting {
	case "msd":
		return m.lessbyMSDs(i, j)
	default:
		return m.lessbyMolIDs(i, j)

	}

}

func (m *mSD) SortBy(sorting string) {
	m.sorting = strings.ToLower(sorting)
	sort.Stable(m)
}

//helper functions

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
func res2atoms(mol chem.Atomer, res []int, names []string, chains []string) []int {
	ret := make([]int, 0, len(res))
	for i := 0; i < mol.Len(); i++ {
		a := mol.Atom(i)

		if (res == nil || isInInt(a.MolID, res)) && isInString(a.Name, names) {
			if chains == nil || isInString(a.Chain, chains) {
				ret = append(ret, i)
			}
		}
	}
	return ret

}

//returns a slice of MSDs where the values with an index not in indexes have been removed.
func msdFilter(msd []float64, indexes []int) []float64 {
	ret := make([]float64, 0, len(indexes))
	for i, v := range msd {
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
