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
	"sort"
	"strings"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/traj/dcd"
	"github.com/rmera/gochem/traj/stf"
	"github.com/rmera/gochem/traj/xtc"
	v3 "github.com/rmera/gochem/v3"
)

//ConcTraje is an interface for a trajectory that can be read concurrently.
type ConcTrajCloser interface {
	chem.ConcTraj
	//closes the trajectory
	Close()
}

type MolIDandChain struct {
	molid int
	chain string
}

func (M *MolIDandChain) String() string {
	return fmt.Sprintf("molid:%04d-chain:%s", M.molid, M.chain)
}

func (M *MolIDandChain) MolID() int {
	return M.molid
}

func (M *MolIDandChain) Chain() string {
	return M.chain
}

//LOVOReturn contains the information returned by LOVOnMostRigid
type LOVOReturn struct {
	N          int
	FramesRead int              //how many frames where actually read from the trajectory
	Natoms     []int            //IDs of the N most rigid atoms
	Nmols      []*MolIDandChain //IDs for the N most rigid residues/molecules
	RMSD       []float64        //the RMSD for all residues, not only the N most rigid.
	Iterations int              //The iterations that were needed for convergency.

}

//String returns a string representation of the LOVOReturn object.
func (L *LOVOReturn) String() string {
	retsl := make([]string, 0, len(L.Nmols)+2)
	retsl = append(retsl, fmt.Sprintf("N: %d,  Frames read: %d, Iterations needed: %d, RMSD: %v, Natoms: %v, Nmols:", L.N, L.FramesRead, L.Iterations, L.RMSD, L.Natoms))
	for _, v := range L.Nmols {
		retsl = append(retsl, v.String())
	}
	return strings.Join(retsl, " ")
}

//Implements Sort to sort by residue
func (L *LOVOReturn) Less(i, j int) bool {
	return L.Nmols[i].molid < L.Nmols[j].molid
}

func (L *LOVOReturn) Swap(i, j int) {
	L.Nmols[i], L.Nmols[j] = L.Nmols[j], L.Nmols[i]
}
func (L *LOVOReturn) Len() int {
	return len(L.Nmols)
}

//PyMOLSel returns a string of text to create a PyMOL
//selection with the L.N most rigid residues from a LOVO
//calculation.
func (L *LOVOReturn) PyMOLSel() string {
	pymolsele := "select rigid,"
	for i, v := range L.Nmols {
		if v.chain == "" {
			pymolsele += fmt.Sprintf(" (resi %d) ", v.molid)

		} else {
			pymolsele += fmt.Sprintf(" (resi %d and chain %s) ", v.molid, v.chain)
		}
		if i < len(L.Nmols)-1 {
			pymolsele += "or"
		}
	}
	return pymolsele

}

//Return the LOVO-selected residues in a format that
//can be pasted into the Gromacs gmx make_ndx tool
//to build a selection. The selection will incluede
//the residues from the chain given, and the chain ID,
//if given, or all residues matching, if not.
//The atomname (normally CA in atomistic simulations or
//BB in Martini ones is not added, mostly out of laziness
//and because it's easy to
//simply append " & a CA" or " & a BB" at the end of the string.
func (L *LOVOReturn) GMX(chain ...string) string {
	data := L.Nmols
	first := true
	ret := ""
	for _, v := range data {
		if len(chain) > 0 && chain[0] != "" && chain[0] != v.Chain() {
			continue
		}
		mid := v.MolID()
		if first {
			ret = fmt.Sprintf("r %d ", mid)
			first = false
			continue
		}
		ret = fmt.Sprintf("%s | r %d ", ret, mid)
	}

	if len(chain) > 0 && chain[0] != "" {
		ret = fmt.Sprintf("%s & chain %s", ret, chain[0])
	}
	return ret
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
	case "stf":
		ret, _, err = stf.New(name)
	case "pdb":
		ret, err = chem.PDBFileRead(name)
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
//If you use this function in your research, please cite the references for the LOVO alignment method:
//10.1371/journal.pone.0119264.
//10.1186/1471-2105-8-306
//If you use this function in a to-consumer program, please print a message asking to cite the references when the function is used.
func LOVO(mol chem.Atomer, ref *v3.Matrix, trajname string, opt ...*Options) (*LOVOReturn, error) {
	var o *Options
	if len(opt) == 0 {
		o = DefaultOptions()
	} else {
		o = opt[0]
	}
	//we first do one iteration aligning the whole thing.
	//these are atom, not residue indexes!
	const defaultTopRigid int = 10 //by default we obtain the 1/topdefaultTopRigid top rigid
	var fullindexes []int = res2atoms(mol, o.atomNames, o.chains, nil)
	var printtraj string
	printtraj = o.writeTraj
	o.writeTraj = ""
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
			o.writeTraj = printtraj
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
		if o.lessThanRMSD > 0 {
			n, newlim = mobileLimit(o.lessThanRMSD, o.minimumN, RMSD)
		}
		indexes := RMSD.MolIDsCopy()
		//There is a difference between the original paper and this implementation.
		//In the original paper, the criterion for convergency is the sum of MSD of the phi*N most rigid atoms
		//Here I simply check that the identity of those phi*N atoms has converged.
		//I think it is completely equivalent.
		if sameElementsInt(indexes[0:n], indexesold[0:n]) && (newlim == o.lessThanRMSD || approxEq(newlim, prevlim, 0.01)) {
			break //converged
		}
		prevlim = newlim
		disagreement = disagreementInt(indexes[0:n], indexesold[0:n])
		indexesold = indexes
	}
	//Now indexes old should countain what we want
	RMSD.SortBy("molid")
	r := &LOVOReturn{N: n, Natoms: indexesold[0:n], Nmols: iDs2MolIDandChains(mol, indexesold[0:n]), RMSD: RMSD.rMSDs, FramesRead: totalframes, Iterations: itercount}
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
	chunksize := o.cpus
	chunk := make([]*v3.Matrix, chunksize)
	tmps := make([]*v3.Matrix, chunksize)
	for i, _ := range chunk {
		chunk[i] = v3.Zeros(mol.Len())
		tmps[i] = v3.Zeros(mol.Len())
	}
	var wtraj *dcd.DCDWObj
	if o.writeTraj != "" {
		wtraj, err = dcd.NewWriter(o.writeTraj, mol.Len())
		if err != nil {
			log.Printf("Couldn't open %s for writing. Will not write aligned trajectory", o.writeTraj) //no error handling here. If we can't open the file, we just don't write anything.
		} else {
			defer wtraj.Close()

		}

	}
	begin := o.begin / chunksize //how many chunks to skip before begining.
	skip := o.skip / chunksize   //how many chunks to skip
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

func isInMIDaC(s *MolIDandChain, c []*MolIDandChain) bool {
	for _, v := range c {
		if v.molid == s.molid && v.chain == s.chain {
			return true
		}
	}
	return false

}

func iDs2MolIDandChains(mol chem.Atomer, indexes []int) []*MolIDandChain {
	ret := make([]*MolIDandChain, 0, len(indexes))
	for _, v := range indexes {
		at := mol.Atom(v)
		s := &MolIDandChain{molid: at.MolID, chain: at.Chain}
		if !isInMIDaC(s, ret) {
			ret = append(ret, s)
		}
	}
	return ret

}
