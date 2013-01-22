/*
 * chem.go, part of gochem.
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

//Package chem provides atom and molecule structures, facilities for reading and writing some
//files used in computational chemistry and some functions for geometric manipulations and shape
//indicators.
package chem

import "fmt"
import "github.com/skelterjohn/go.matrix"

//import "strings"

/**Note: Many funcitons here panic instead of returning errors. This is because they are "fundamental" 
 * functions. I considered that if something goes wrong here, the program is way-most likely wrong and should
 * crash. Most panics are related to using the funciton on a nil object or trying to access out-of bounds 
 * fields**/

//Atom contains the atoms read except for the coordinates, which will be in a matrix
//and the b-factors, which are in a separate slice of float64.
type Atom struct {
	Name      string
	Id        int
	Tag       int //Just added this for something that someone might want to keep that is not a float.
	Molname   string
	Molname1  byte //the one letter name for residues and nucleotids
	Molid     int
	Chain     byte
	Mass      float64 //hopefully all these float64 are not too much memory
	Occupancy float64
	Vdw       float64
	Charge    float64
	Symbol    string
	Het       bool // is hetatm in the pdb file?
}

//Atom methods

//Copy returns a copy of the Atom object.
func (A *Atom) Copy() *Atom {
	if A == nil {
		panic("Attempted to copy a nil atom")
	}
	Newat := new(Atom)
	Newat.Name = A.Name
	Newat.Id = A.Id
	Newat.Tag = A.Tag
	Newat.Molname = A.Molname
	Newat.Molname1 = A.Molname1
	Newat.Molid = A.Molid
	Newat.Chain = A.Chain
	Newat.Mass = A.Mass
	Newat.Occupancy = A.Occupancy
	Newat.Vdw = A.Vdw
	Newat.Charge = A.Charge
	Newat.Symbol = A.Symbol
	Newat.Het = A.Het
	return Newat
}

/*****Topology type***/

//Topology contains information about a molecule which is not expected to change in time (i.e. everything except for coordinates and b-factors)
type Topology struct {
	Atoms    []*Atom
	charge   int
	unpaired int
}

//MakeTopology makes a molecule with ats atoms, coords coordinates, bfactors b-factors 
//charge charge and unpaired unpaired electrons, and returns it. It returns error if 
//one of the slices is nil. It doesnt check for consitensy across slices or correct charge
//or unpaired electrons.
func MakeTopology(ats []*Atom, charge, unpaired int) (*Topology, error) {
	if ats == nil {
		return nil, fmt.Errorf("Supplied a nil Topology")
	}
	top := new(Topology)
	top.Atoms = ats
	top.charge = charge
	top.unpaired = unpaired
	return top, nil
}

/*Topology methods*/

//Charge gets the total charge of the topology
func (T *Topology) Charge() int {
	return T.charge
}

//Unpaired gets the number of unpaired electrons in the topology
func (T *Topology) Unpaired() int {
	return T.unpaired
}

//SetCharge sets the total charge of the topology to i
func (T *Topology) SetCharge(i int) {
	T.charge = i
}

//SetUnpaired sets the number of unpaired electrons in the topology to i
func (T *Topology) SetUnpaired(i int) {
	T.unpaired = i
}



//Sets the current order of atoms as Id and the order of molecules as
//Molid for all atoms
func (T *Topology) ResetIds(){
	currid:=1
	currid2:=1
	for key,val:=range(T.Atoms){
		T.Atoms[key].Id=key+1
		if currid==val.Molid{
			continue
			}
		if currid==val.Molid-1{ //We hit a new molecule
			currid2++
			currid++
			continue
			}
		//change of residue after fixing one that didnt match position
		if currid2!=val.Molid{
			currid2=T.Atoms[key].Molid
			T.Atoms[key].Molid=currid+1
			currid=currid+1
			continue
		}
		//A residue's index doesnt match its position
		T.Atoms[key].Molid=currid
		
		
		
	}
}


//Copy returns a copy of the molecule including coordinates
func (T *Topology) CopyAtoms() Ref {
	Top := new(Topology)
	Top.Atoms = make([]*Atom, T.Len())
	for key, val := range T.Atoms {
		Top.Atoms[key] = val.Copy()
	}
	Top.charge = T.charge
	Top.unpaired = T.unpaired
	return Top
}

//Atom returns the Atom corresponding to the index i
//of the Atom slice in the Topology. Panics if 
//out of range.
func (T *Topology) Atom(i int) *Atom {
	if i >= T.Len() {
		panic("Topology: Requested Atom out of bounds")
	}
	return T.Atoms[i]
}

//SetAtom sets the (i+1)th Atom of the topology to aT.
//Panics if out of range
func (T *Topology) SetAtom(i int, at *Atom) {
	if i >= T.Len() {
		panic("Topology: Tried to set Atom out of bounds")
	}
	T.Atoms[i] = at
}

//AddAtom appends an atom at the end of the reference
func (T *Topology) AddAtom(at *Atom) Ref {
	newmol, ok := T.CopyAtoms().(*Topology)
	if ok == false {
		panic("cant happen")
	}
	newmol.Atoms = append(newmol.Atoms, at)
	return newmol
}

//SelectAtoms, given a list of ints,  returns an array of the atoms with the
//corresponding position in the molecule
//Changes to these atoms affect the original reference.
//The charge and multiplicity (unpaired electrons) for the molecule is just the one
//for the parent reference and its not guarranteed to be correct.
func (T *Topology) SomeAtoms(atomlist []int) (Ref, error) {
	var err error
	var ret []*Atom
	lenatoms := len(T.Atoms)
	for k, j := range atomlist {
		if j > lenatoms-1 {
			return nil, fmt.Errorf("Atom requested (Number: %d, value: %d) out of range!", k, j)
		}
		ret = append(ret, T.Atoms[j])
	}
	finalret,err:=MakeTopology(ret, T.Charge(), T.Unpaired())
	return finalret, err
}

//Returns a copy of T with the i atom deleted by reslicing
//this means that the copy still uses as much memory as the original T!
func (T *Topology) DelAtom(i int) Ref {
	if i >= T.Len() {
		panic("Topology: Tried to delete Atom out of bounds")
	}
	N, ok := T.CopyAtoms().(*Topology)
	if ok == false {
		panic("cant happen!")
	}
	if i == N.Len()-1 {
		N.Atoms = N.Atoms[:i]
	} else {
		N.Atoms = append(N.Atoms[:i], N.Atoms[i+1:]...)
	}
	return N
}

//Len returns the length of the molecule.
func (T *Topology) Len() int {
	return len(T.Atoms) //This shouldnt be needed
}

//MassCol returns a DenseMatrix 1-col matrix with masses of atoms and an error if they have not been calculated
func (T *Topology) MassCol() (*matrix.DenseMatrix, error) {
	mass := make([]float64, T.Len())
	for i := 0; i < T.Len(); i++ {
		thisatom := T.Atom(i)
		if thisatom.Mass == 0 {
			return nil, fmt.Errorf("Not all the masses have been obtained: %d %v", i, thisatom)
		}
		mass[i] = thisatom.Mass
	}
	massmat := matrix.MakeDenseMatrix(mass, len(mass), 1)
	return massmat, nil
}

/**Type Molecule**/

//Molecule contains all the info for a molecule in many states. The info that is expected to change between states,
//Coordinates and b-factors are stored separately from other atomic info.
type Molecule struct {
	*Topology
	Coords   []*matrix.DenseMatrix
	Bfactors []*matrix.DenseMatrix
	current  int
}

//MakeMolecule makes a molecule with ats atoms, coords coordinates, bfactors b-factors 
//charge charge and unpaired unpaired electrons, and returns it. It returns error if 
//one of the slices is nil. It doesnt check for consitensy across slices or correct charge
//or unpaired electrons.
func MakeMolecule(ats Ref, coords, bfactors []*matrix.DenseMatrix) (*Molecule, error) {
	if ats == nil {
		return nil, fmt.Errorf("Supplied a nil Reference")
	}
	if coords == nil {
		return nil, fmt.Errorf("Supplied a nil Coords slice")
	}
	if bfactors == nil {
		return nil, fmt.Errorf("Supplied a nil Bfactors slice")
	}
	mol := new(Molecule)
	top, ok := ats.(*Topology)
	if ok == true {
		mol.Topology = top
	} else {
		mol.Topology = new(Topology)
		mol.Atoms = make([]*Atom, ats.Len())
		for i := 0; i < ats.Len(); i++ {
			mol.Atoms[i] = ats.Atom(i)
		}
	}
	mol.Coords = coords
	mol.Bfactors = bfactors
	return mol, nil
}

//The molecule methods:

//Deletes the coodinate i from every frame of the molecule.
func (M *Molecule) DelCoord(i int) error {
	var err error
	if i >= M.Coords[0].Rows() {
		panic("Reference: Tried to delete Coord out of bounds")
	}
	for j := 0; j < len(M.Coords); j++ {
		M.Coords[j], err = DMDelRow(M.Coords[j], i)
		if err != nil {
			return err
		}
	}
	return nil
}

//Deletes atom i and its coordinates from the molecule.
func (M *Molecule) Del(i int) error {
	if i >= M.Len() {
		panic("Tried to delete Atom out of bounds")
	}
	if i == M.Len()-1 {
		M.Atoms = M.Atoms[:i]
	} else {
		M.Atoms = append(M.Atoms[:i], M.Atoms[i+1:]...)
	}
	err := M.DelCoord(i)
	return err
}

//Copy returns a copy of the molecule including coordinates
func (M *Molecule) Copy() *Molecule {
	if err := M.Corrupted(); err != nil {
		panic(err.Error())
	}
	mol := new(Molecule)
	/*The foll. code is from Topology.CopyAtoms, I wrote it again just not to deal with the type assersion 
	since CopyAtoms returns a Ref interface. This should change at some point*/
	mol.Topology = new(Topology)
	mol.Atoms = make([]*Atom, M.Len())
	for key, val := range M.Atoms {
		mol.Atoms[key] = val.Copy()
	}
	mol.charge = M.charge
	mol.unpaired = M.unpaired
	//Up to previous line
	mol.Coords = make([]*matrix.DenseMatrix, 0, len(M.Coords))
	mol.Bfactors = make([]*matrix.DenseMatrix, 0, len(M.Bfactors))
	for key, val := range M.Coords {
		mol.Coords = append(mol.Coords, val.Copy())
		mol.Bfactors = append(mol.Bfactors, M.Bfactors[key].Copy())
	}
	if err := mol.Corrupted(); err != nil {
		panic(fmt.Sprintf("Molecule creation error: %s", err.Error())) //copying a corrupted molecule means that the program is wrong.
	}
	return mol
}

//AddFrame akes a matrix of coordinates and appends them at the end of the Coords.
// It checks that the number of coordinates matches the number of atoms.
func (M *Molecule) AddFrame(newframe *matrix.DenseMatrix) {
	if newframe == nil {
		panic("Attempted to add nil frame")
	}
	if newframe.Cols() != 3 {
		panic("Malformed coord matrix!")
	}
	if M.Len() != newframe.Rows() {
		panic(fmt.Sprintf("Wrong number of coordinates (%d)", newframe.Rows()))
	}
	if M.Coords == nil {
		M.Coords = make([]*matrix.DenseMatrix, 1, 1)
	}
	M.Coords = append(M.Coords, newframe)
}

//AddManyFrames adds the array of matrices newfames to the molecule. It checks that
//the number of coordinates matches the number of atoms.
func (M *Molecule) AddManyFrames(newframes []*matrix.DenseMatrix) {
	if newframes == nil {
		panic("Attempted to add nil frames")
	}
	atomslen := M.Len()
	if M.Coords == nil {
		M.Coords = make([]*matrix.DenseMatrix, 1, len(newframes))
	}
	//Must add something here to change
	for i := range newframes {
		if atomslen != newframes[i].Rows() {
			panic(fmt.Sprintf("Wrong number of coordinates (%d)  in frame %d", newframes[i].Rows(), i))
		}
		M.Coords = append(M.Coords, newframes[i]) //At the end there shouldnt be so much copying since
		//only pointers are copied.
	}
}

//Coord returns the coords for the atom atom in the frame frame.
//panics if frame or coords are out of range.
func (M *Molecule) Coord(atom, frame int) *matrix.DenseMatrix {
	if frame >= len(M.Coords) {
		panic(fmt.Sprintf("Frame requested (%d) out of range", frame))
	}
	if atom >= M.Coords[frame].Rows() {
		panic(fmt.Sprintf("Requested coordinate (%d) out of bounds (%d)", atom, M.Coords[frame].Rows()))
	}
	return M.Coords[frame].GetRowVector(atom)
}

//Current returns the number of the next readed frame
func (M *Molecule) Current() int {
	if M == nil {
		return -1
	}
	return M.current
}

//SetCurrent sets the value of the frame nex to be read
//to i.
func (M *Molecule) SetCurrent(i int) {
	if i < 0 || i >= len(M.Coords) {
		panic("Invalid new value for current")
	}
	M.current = i
}

//SetCoords replaces the coordinates of atoms in the positions given by atomlist with the ones in newcoords (in order)
//If atomlist contains a single element, it replaces as many coordinates as given in newcoords, starting 
//at the element in atomlist. In the latter case, the function checks that there are enough coordinates to
//replace and returns an error if not.
func (M *Molecule) SetCoords(newcoords *matrix.DenseMatrix, atomlist []int, frame int) {
	if frame >= len(M.Coords) {
		panic(fmt.Sprintf("Frame (%d) out of range!", frame))
	}
	//If supplies a list with one number, the newcoords will replace the old coords
	//Starting that number. We do check that you don't put more coords than spaces we have.
	if len(atomlist) == 1 {
		if newcoords.Rows() > M.Coords[frame].Rows()-atomlist[0]-1 {
			panic(fmt.Sprintf("Cant replace starting from position %d: Not enough atoms in molecule", atomlist[0]))
		}
		M.Coords[frame].SetMatrix(atomlist[0], 0, newcoords)
		return
	}
	//If the list has more than one atom
	lenatoms := M.Len()
	for k, j := range atomlist {
		if j > lenatoms-1 {
			panic(fmt.Sprintf("Requested position number: %d (%d) out of range", k, j))
		}
		M.Coords[frame].SetMatrix(j, 0, newcoords.GetRowVector(k))
	}
}

//Corrupted checks whether the molecule is corrupted, i.e. the
//coordinates don't match the number of atoms. It also checks
//That the coordinate matrices have 3 columns.
func (M *Molecule) Corrupted() error {
	var err error
	if M.Bfactors == nil {
		M.Bfactors = make([]*matrix.DenseMatrix, 0, len(M.Coords))
		M.Bfactors = append(M.Bfactors, matrix.Zeros(M.Len(), 1))
	}
	lastbfac := len(M.Bfactors) - 1
	for i := range M.Coords {
		if M.Len() != M.Coords[i].Rows() || M.Coords[i].Cols() != 3 {
			err = fmt.Errorf("Inconsistent coordinates/atoms in frame %d: Atoms %d, coords: %d", i, M.Len(), M.Coords[i].Rows())
			break
		}
		//Since bfactors are not as important as coordinates, we will just fill with 
		//zeroes anything that is lacking or incomplete instead of returning an error.
		if lastbfac < i {
			bfacs := matrix.Zeros(M.Len(), 1)
			M.Bfactors = append(M.Bfactors, bfacs)
		} else if M.Bfactors[i].Rows() < M.Len() {
			M.Bfactors[i] = matrix.Zeros(M.Len(), 1)
		}
	}
	return err
}

//LenFrames returns the number of frames in the molecule
func (M *Molecule) LenFrames() int {
	return len(M.Coords)
}

//Implementaiton of the sort.Interface

//Swap function, as demanded by sort.Interface. It swaps atoms, coordinates 
//(all frames) and bfactors of the molecule.	
func (M *Molecule) Swap(i, j int) {
	M.Atoms[i], M.Atoms[j] = M.Atoms[j], M.Atoms[i]
	for k := 0; k < len(M.Coords); k++ {
		M.Coords[k].SwapRows(i, j)
		M.Bfactors[k].SwapRows(i, j)
	}
}

//Less: Should the atom i be sorted before atom j?
func (M *Molecule) Less(i, j int) bool {
	return M.Bfactors[0].Get(i, 0) < M.Bfactors[0].Get(j, 0)
}

//Len is implemented in Topology

//End sort.Interface

/******************************************
//The following implement the Traj interface
**********************************************/

//Checks that the molecule exists and has some existent
//Coordinates, in which case returns true.
//It returns false otherwise.
// The coordinates could still be empty
func (M *Molecule) Readable() bool {
	if M != nil || M.Coords != nil ||  M.current < len(M.Coords){
		
		return true
	}
	return false
}

//Returns the  next frame and an error
func (M *Molecule) Next(a bool) (*matrix.DenseMatrix, error) {
	if M.current >= len(M.Coords) {
		return nil, fmt.Errorf("No more frames")
	}
	if a == false {
		return nil, nil
	}
	M.current++
	return M.Coords[M.current-1], nil
}

//Initializes molecule to be read as a traj (not tested!)
func (M *Molecule) InitRead() error {
	if M==nil || len(M.Coords)==0 {
		return fmt.Errorf("Bad molecule")	
	}
	M.current=0
	return nil
}

/*NextConc takes a slice of bools and reads as many frames as elements the list has
form the trajectory. The frames are discarted if the corresponding elemetn of the slice
* is false. The function returns a slice of channels through each of each of which 
* a *matrix.DenseMatrix will be transmited*/
func (M *Molecule) NextConc(frames []bool) ([]chan *matrix.DenseMatrix, error) {
	toreturn := make([]chan *matrix.DenseMatrix, 0, len(frames))
	used := false
	for _, val := range frames {
		if val == false {
			M.current++
			toreturn = append(toreturn, nil)
			continue
		}
		if M.current >= len(M.Coords) {
			if used == false {
				return nil, fmt.Errorf("No more frames")
			} else {
				return toreturn, fmt.Errorf("No more frames")
			}
		}
		used = true
		toreturn = append(toreturn, make(chan *matrix.DenseMatrix))
		go func(a *matrix.DenseMatrix, pipe chan *matrix.DenseMatrix) {
			pipe <- a
		}(M.Coords[M.current], toreturn[len(toreturn)-1])
		M.current++
	}
	return toreturn, nil
}



/**End Traj interface implementation***********/

//End Molecule methods
