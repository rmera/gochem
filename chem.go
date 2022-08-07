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

package chem

import (
	"fmt"
	"sort"

	v3 "github.com/rmera/gochem/v3"
)

//import "strings"

/* Many funcitons here panic instead of returning errors. This is because they are "fundamental"
 * functions. I considered that if something goes wrong here, the program is way-most likely wrong and should
 * crash. Most panics are related to using the funciton on a nil object or trying to access out-of bounds
 * fields
 */

//Atom contains the information to represent an atom, except for the coordinates, which will be in a separate *v3.Matrix
//and the b-factors, which are in a separate slice of float64.
type Atom struct {
	Name      string  //PDB name of the atom
	ID        int     //The PDB index of the atom
	index     int     //The place of the atom in a set. I won't make it accessible to ensure that it does correspond to the ordering.
	Tag       int     //Just added this for something that someone might want to keep that is not a float.
	MolName   string  //PDB name of the residue or molecule (3-letter code for residues)
	MolName1  byte    //the one letter name for residues and nucleotids
	Char16    byte    //Whatever is in the column 16 (counting from 0) in a PDB file, anything.
	MolID     int     //PDB index of the corresponding residue or molecule
	Chain     string  //One-character PDB name for a chain.
	Mass      float64 //hopefully all these float64 are not too much memory
	Occupancy float64 //a PDB crystallographic field, often used to store values of interest.
	Vdw       float64 //radius
	Charge    float64 //Partial charge on an atom
	Symbol    string
	Het       bool    // is the atom an hetatm in the pdb file? (if applicable)
	Bonds     []*Bond //The bonds connecting the atom to others.
}

//Atom methods

//Copy returns a copy of the Atom object.
//puts the copy into the
func (N *Atom) Copy(A *Atom) {
	if A == nil || N == nil {
		panic(ErrNilAtom)
	}
	N.Name = A.Name
	N.ID = A.ID
	N.Tag = A.Tag
	N.MolName = A.MolName
	N.MolName1 = A.MolName1
	N.MolID = A.MolID
	N.Chain = A.Chain
	N.Mass = A.Mass
	N.Occupancy = A.Occupancy
	N.Vdw = A.Vdw
	N.Charge = A.Charge
	N.Symbol = A.Symbol
	N.Het = A.Het
}

//Index returns the index of the atom
func (N *Atom) Index() int {
	return N.index
}

/*****Topology type***/

//Topology contains information about a molecule which is not expected to change in time (i.e. everything except for coordinates and b-factors)
type Topology struct {
	Atoms  []*Atom
	charge int
	multi  int
}

//NewTopology returns topology with ats atoms,
//charge charge and multi multiplicity.
// It doesnt check for consitency across slices, correct charge
//or unpaired electrons.
func NewTopology(charge, multi int, ats ...[]*Atom) *Topology {
	top := new(Topology)
	if len(ats) == 0 || ats[0] == nil {
		top.Atoms = make([]*Atom, 0, 0) //return nil, fmt.Errorf("Supplied a nil Topology")
	} else {
		top.Atoms = ats[0]
	}
	top.charge = charge
	top.multi = multi
	return top
}

/*Topology methods*/

//Charge returns the total charge of the topology
func (T *Topology) Charge() int {
	return T.charge
}

//Multi returns the multiplicity in the topology
func (T *Topology) Multi() int {
	return T.multi
}

//SetCharge sets the total charge of the topology to i
func (T *Topology) SetCharge(i int) {
	T.charge = i
}

//SetMulti sets the multiplicity in the topology to i
func (T *Topology) SetMulti(i int) {
	T.multi = i
}

//FillMasses tries to get fill the  masses for atom that don't have one
//by getting it from the symbol. Only a few common elements are supported
func (T *Topology) FillMasses() {
	for _, val := range T.Atoms {
		if val.Symbol != "" && val.Mass == 0 {
			val.Mass = symbolMass[val.Symbol] //Not error checking
		}
	}
}

//FillsIndexes sets the Index value of each atom to that cooresponding to its
//place in the molecule.
func (T *Topology) FillIndexes() {
	for key, val := range T.Atoms {
		val.index = key
	}

}

//FillVdw tries to get fill the  van der Waals radii for the atoms in the molecule
//from a symbol->radii map. Only a few common elements are supported
func (T *Topology) FillVdw() {
	for _, val := range T.Atoms {
		if val.Symbol != "" && val.Vdw == 0 {
			val.Vdw = symbolVdwrad[val.Symbol] //Not error checking
		}
	}
}

//ResetIDs sets the current order of atoms as ID and the order of molecules as
//MolID for all atoms
func (T *Topology) ResetIDs() {
	currid := 1
	currid2 := 1
	for key, val := range T.Atoms {
		T.Atoms[key].ID = key + 1
		if currid == val.MolID {
			continue
		}
		if currid == val.MolID-1 { //We hit a new molecule
			currid2++
			currid++
			continue
		}
		//change of residue after fixing one that didnt match position
		if currid2 != val.MolID {
			currid2 = T.Atoms[key].MolID
			T.Atoms[key].MolID = currid + 1
			currid = currid + 1
			continue
		}
		//A residue's index doesnt match its position
		T.Atoms[key].MolID = currid

	}
}

//CopyAtoms copies the atoms form A into the receiver topology. This is a deep copy, so the receiver must have
//at least as many atoms as A.
func (T *Topology) CopyAtoms(A Atomer) {
	//T := new(Topology)
	T.Atoms = make([]*Atom, A.Len())
	for key := 0; key < A.Len(); key++ {
		T.Atoms[key].Copy(A.Atom(key))
	}
}

//Atom returns the Atom corresponding to the index i
//of the Atom slice in the Topology. Panics if
//out of range.
func (T *Topology) Atom(i int) *Atom {
	if i >= T.Len() {
		panic(ErrAtomOutOfRange)
	}
	return T.Atoms[i]
}

//SetAtom sets the (i+1)th Atom of the topology to aT.
//Panics if out of range
func (T *Topology) SetAtom(i int, at *Atom) {
	if i >= T.Len() {
		panic(ErrAtomOutOfRange)
	}
	T.Atoms[i] = at
}

//AppendAtom appends an atom at the end of the reference
func (T *Topology) AppendAtom(at *Atom) {
	/*newmol, ok := T.CopyAtoms().(*Topology)
	if !ok {
		panic("cant happen")
	}
	newmol.Atoms = append(newmol.Atoms, at)*/
	T.Atoms = append(T.Atoms, at)
}

//SelectAtoms puts the subset of atoms in T that have
//indexes in atomlist into the receiver. Panics if problem.
func (R *Topology) SomeAtoms(T Atomer, atomlist []int) {
	var ret []*Atom
	lenatoms := T.Len()
	for k, j := range atomlist {
		if j > lenatoms-1 {
			err := fmt.Sprintf("goChem: Atom requested (Number: %d, value: %d) out of range", k, j)
			panic(PanicMsg(err))
		}
		ret = append(ret, T.Atom(j))
	}
	R.Atoms = ret
}

//SelectAtomsSafe puts the atoms of T
//with indexes in atomlist into the receiver. Returns error if problem.
func (R *Topology) SomeAtomsSafe(T Atomer, atomlist []int) error {
	f := func() { R.SomeAtoms(T, atomlist) }
	return gnMaybe(gnPanicker(f))
}

//DelAtom Deletes atom i by reslicing.
//This means that the copy still uses as much memory as the original T.
func (T *Topology) DelAtom(i int) {
	if i >= T.Len() {
		panic(ErrAtomOutOfRange)
	}
	if i == T.Len()-1 {
		T.Atoms = T.Atoms[:i]
	} else {
		T.Atoms = append(T.Atoms[:i], T.Atoms[i+1:]...)
	}
}

//Len returns the number of atoms in the topology.
func (T *Topology) Len() int {
	//if T.Atoms is nil, return len(T.Atoms) will panic, so I will let that happen for now.
	//	if T.Atoms == nil {
	//		panic(ErrNilAtoms)
	//	}
	return len(T.Atoms)
}

//MassCol returns a slice of float64 with the masses of the atoms in the topology, or nil and an error if they have not been calculated
func (T *Topology) Masses() ([]float64, error) {
	mass := make([]float64, T.Len())
	for i := 0; i < T.Len(); i++ {
		thisatom := T.Atom(i)
		if thisatom.Mass == 0 {
			return nil, CError{fmt.Sprintf("goChem: Not all the masses have been obtained: %d %v", i, thisatom), []string{"Topology.Masses"}}
		}
		mass[i] = thisatom.Mass
	}
	return mass, nil
}

//AssignBonds assigns bonds to a molecule based on a simple distance
//criterium, similar to that described in DOI:10.1186/1758-2946-3-33
func (T *Topology) AssignBonds(coord *v3.Matrix) error {
	// might get slow for
	//large systems. It's really not thought
	//for proteins or macromolecules.
	//For this reason, this method is not called automatically when building a new topology.
	//Well, that and that it requires a Matrix object, which would mean changing the
	//signature of the NewTopology function.
	var t1, t2 *v3.Matrix
	var at1, at2 *Atom
	T.FillIndexes()
	t3 := v3.Zeros(1)
	bonds := make([]*Bond, 0, 10)
	tot := T.Len()
	var nextIndex int
	for i := 0; i < tot; i++ {
		t1 = coord.VecView(i)
		at1 = T.Atoms[i]
		cov1 := symbolCovrad[at1.Symbol]
		if cov1 == 0 {
			err := new(CError)
			err.msg = fmt.Sprintf("Couldn't find the covalent radii  for %s %d", at1.Symbol, i)
			err.Decorate("AssignBonds")
			return err
		}
		for j := i + 1; j < tot; j++ {
			t2 = coord.VecView(j)
			at2 = T.Atoms[j]
			cov2 := symbolCovrad[at2.Symbol]
			if cov2 == 0 {
				err := new(CError)
				err.msg = fmt.Sprintf("Couldn't find the covalent radii  for %s %d", at2.Symbol, j)
				err.Decorate("AssignBonds")
				return err
			}

			t3.Sub(t2, t1)
			d := t3.Norm(2)
			if d < cov1+cov2+bondtol && d > tooclose {
				b := &Bond{Index: nextIndex, Dist: d, At1: at1, At2: at2}
				at1.Bonds = append(at1.Bonds, b)
				at2.Bonds = append(at2.Bonds, b)
				bonds = append(bonds, b) //just to easily keep track of them.
				nextIndex++
			}

		}
	}

	//Now we check that no atom has too many bonds.
	for i := 0; i < tot; i++ {
		at := T.Atoms[i]
		max := symbolMaxBonds[at.Symbol]
		if max == 0 { //means there is not a specified number of bonds for this atom.
			continue
		}
		sort.Slice(at.Bonds, func(i, j int) bool { return at.Bonds[i].Dist < at.Bonds[j].Dist })
		//I am hoping this will remove bonds until len(at.Bonds) is not
		//greater than max.
		for i := len(at.Bonds); i > max; i = len(at.Bonds) {
			err := at.Bonds[len(at.Bonds)-1].Remove() //we remove the longest bond
			if err != nil {
				return errDecorate(err, "AssignBonds")
			}
		}

	}

	return nil
}

/**Type Molecule**/

//Molecule contains all the info for a molecule in many states. The info that is expected to change between states,
//Coordinates and b-factors are stored separately from other atomic info.
type Molecule struct {
	*Topology
	Coords      []*v3.Matrix
	Bfactors    [][]float64
	XYZFileData []string //This can be anything. The main rationale for including it is that XYZ files have a "comment"
	//line after the first one. This line is sometimes used to write the energy of the structure.
	//So here the line can be kept for each XYZ frame, and parse later
	current int
}

//NewMolecule makes a molecule with ats atoms, coords coordinates, bfactors b-factors
//charge charge and unpaired unpaired electrons, and returns it. It doesnt check for
// consitency across slices or correct charge or unpaired electrons.
func NewMolecule(coords []*v3.Matrix, ats Atomer, bfactors [][]float64) (*Molecule, error) {
	if ats == nil {
		return nil, CError{"Supplied a nil Reference", []string{"NewMolCule"}}
	}
	if coords == nil {
		return nil, CError{"Supplied a nil Coord slice", []string{"NewMolCule"}}
	}
	//	if bfactors == nil {
	//		return nil, fmt.Errorf("Supplied a nil Bfactors slice")
	//	}
	mol := new(Molecule)
	atcopy := func() {
		mol.Topology = NewTopology(9999, -1, make([]*Atom, 0, ats.Len())) //I use 9999 for charge and -1 or multi to indicate that they are not truly set. So far NewTopology never actually returns any error so it's safe to ignore them.
		for i := 0; i < ats.Len(); i++ {
			mol.Atoms = append(mol.Atoms, ats.Atom(i))
		}
	}
	switch ats := ats.(type) { //for speed
	case *Topology:
		mol.Topology = ats
	case AtomMultiCharger:
		atcopy()
		mol.SetMulti(ats.Multi())
		mol.SetCharge(ats.Charge())
	default:
		atcopy()
	}
	mol.Coords = coords
	mol.Bfactors = bfactors
	return mol, nil
}

//The molecule methods:

//DelCoord deletes the coodinate i from every frame of the molecule.
func (M *Molecule) DelCoord(i int) error {
	//note: Maybe this shouldn't be exported. Unexporting it could be a reasonable API change.
	r, _ := M.Coords[0].Dims()
	var err error
	for j := 0; j < len(M.Coords); j++ {
		tmp := v3.Zeros(r - 1)
		tmp.DelVec(M.Coords[j], i)
		M.Coords[j] = tmp
		if err != nil {
			return err
		}
	}
	return nil
}

//Del Deletes atom i and its coordinates from the molecule.
func (M *Molecule) Del(i int) error {
	if i >= M.Len() {
		panic(ErrAtomOutOfRange)
	}
	M.DelAtom(i)
	err := M.DelCoord(i)
	return err
}

//Copy puts in the receiver a copy of the molecule  A including coordinates
func (M *Molecule) Copy(A *Molecule) {
	if err := A.Corrupted(); err != nil {
		panic(err.Error())
	}
	r, _ := A.Coords[0].Dims()
	//mol := new(Molecule)
	M.Topology = new(Topology)
	for i := 0; i < A.Len(); i++ {
		at := new(Atom)
		at.Copy(A.Atom(i))
		M.Topology.Atoms = append(M.Topology.Atoms, at)

	}
	//	M.CopyAtoms(A)
	M.Coords = make([]*v3.Matrix, 0, len(A.Coords))
	M.Bfactors = make([][]float64, 0, len(A.Bfactors))
	for key, val := range A.Coords {
		tmp := v3.Zeros(r)
		tmp.Copy(val)
		M.Coords = append(M.Coords, tmp)
		tmp2 := copyB(A.Bfactors[key])
		M.Bfactors = append(M.Bfactors, tmp2)
	}
	if err := M.Corrupted(); err != nil {
		panic(PanicMsg(fmt.Sprintf("goChem: Molecule creation error: %s", err.Error())))
	}
}

func copyB(b []float64) []float64 {
	r := make([]float64, len(b), len(b))
	for k, v := range b {
		r[k] = v
	}
	return r
}

//AddFrame akes a matrix of coordinates and appends them at the end of the Coords.
// It checks that the number of coordinates matches the number of atoms.
func (M *Molecule) AddFrame(newframe *v3.Matrix) {
	if newframe == nil {
		panic(ErrNilFrame)
	}
	r, c := newframe.Dims()
	if c != 3 {
		panic(ErrNotXx3Matrix)
	}
	if M.Len() != r {
		panic(PanicMsg(fmt.Sprintf("goChem: Wrong number of coordinates (%d)", newframe.NVecs())))
	}
	if M.Coords == nil {
		M.Coords = make([]*v3.Matrix, 1, 1)
	}
	M.Coords = append(M.Coords, newframe)
}

//AddManyFrames adds the array of matrices newfames to the molecule. It checks that
//the number of coordinates matches the number of atoms.
func (M *Molecule) AddManyFrames(newframes []*v3.Matrix) {
	if newframes == nil {
		panic(ErrNilFrame)
	}
	if M.Coords == nil {
		M.Coords = make([]*v3.Matrix, 1, len(newframes))
	}
	for key, val := range newframes {
		f := func() { M.AddFrame(val) }
		err := gnMaybe(gnPanicker(f))
		if err != nil {
			panic(PanicMsg(fmt.Sprintf("goChem: %s in frame %d", err.Error(), key)))
		}
	}
}

//Coord returns the coords for the atom atom in the frame frame.
//panics if frame or coords are out of range.
func (M *Molecule) Coord(atom, frame int) *v3.Matrix {
	if frame >= len(M.Coords) {
		panic(PanicMsg(fmt.Sprintf("goChem: Frame requested (%d) out of range", frame)))
	}
	r, _ := M.Coords[frame].Dims()
	if atom >= r {
		panic(PanicMsg(fmt.Sprintf("goChem: Requested coordinate (%d) out of bounds (%d)", atom, M.Coords[frame].NVecs())))
	}
	ret := v3.Zeros(1)
	empt := M.Coords[frame].VecView(atom)
	ret.Copy(empt)
	return ret
}

//Current returns the number of the next read frame
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
		panic(PanicMsg(fmt.Sprintf("goChem: Invalid new value for current frame: %d Current frames: %d", i, len(M.Coords))))
	}
	M.current = i
}

/*
//SetCoords replaces the coordinates of atoms in the positions given by atomlist with the gnOnes in NewVecs (in order)
//If atomlist contains a single element, it replaces as many coordinates as given in NewVecs, starting
//at the element in atomlist. In the latter case, the function checks that there are enough coordinates to
//replace and returns an error if not.
func (M *Molecule) SetCoords(NewVecs *CoordMA, atomlist []int, frame int) {
	if frame >= len(M.Coords) {
		panic(fmt.Sprintf("Frame (%d) out of range!", frame))
	}
	//If supplies a list with one number, the NewVecs will replace the old coords
	//Starting that number. We do check that you don't put more coords than spaces we have.
	if len(atomlist) == 1 {
		if NewVecs.Rows() > M.Coords[frame].Rows()-atomlist[0]-1 {
			panic(fmt.Sprintf("Cant replace starting from position %d: Not enough atoms in molecule", atomlist[0]))
		}
		M.Coords[frame].SetMatrix(atomlist[0], 0, NewVecs)
		return
	}
	//If the list has more than one atom
	lenatoms := M.Len()
	for k, j := range atomlist {
		if j > lenatoms-1 {
			panic(fmt.Sprintf("Requested position number: %d (%d) out of range", k, j))
		}
		M.Coords[frame].SetMatrix(j, 0, NewVecs.GetRowVector(k))
	}
}

*/

//Corrupted checks whether the molecule is corrupted, i.e. the
//coordinates don't match the number of atoms. It also checks
//That the coordinate matrices have 3 columns.
func (M *Molecule) Corrupted() error {
	var err error
	if M.Bfactors == nil {
		M.Bfactors = make([][]float64, 0, len(M.Coords))
		M.Bfactors = append(M.Bfactors, make([]float64, M.Len()))
	}
	lastbfac := len(M.Bfactors) - 1
	for i := range M.Coords {
		r, c := M.Coords[i].Dims()
		if M.Len() != r || c != 3 {
			err = CError{fmt.Sprintf("Inconsistent coordinates/atoms in frame %d: Atoms %d, coords: %d", i, M.Len(), M.Coords[i].NVecs()), []string{"Molecule.Corrupted"}}
			break
		}
		//Since bfactors are not as important as coordinates, we will just fill with
		//zeroes anything that is lacking or incomplete instead of returning an error.

		if lastbfac < i {
			bfacs := make([]float64, M.Len())
			M.Bfactors = append(M.Bfactors, bfacs)
		}
		bfr := len(M.Bfactors[i])
		if bfr < M.Len() {
			M.Bfactors[i] = make([]float64, M.Len(), 1)
		}
	}
	return err
}

//NFrames returns the number of frames in the molecule
func (M *Molecule) NFrames() int {
	return len(M.Coords)
}

//Implementaiton of the sort.Interface

//Swap function, as demanded by sort.Interface. It swaps atoms, coordinates
//(all frames) and bfactors of the molecule.
func (M *Molecule) Swap(i, j int) {
	M.Atoms[i], M.Atoms[j] = M.Atoms[j], M.Atoms[i]
	for k := 0; k < len(M.Coords); k++ {
		M.Coords[k].SwapVecs(i, j)
		t1 := M.Bfactors[k][i]
		t2 := M.Bfactors[k][j]
		M.Bfactors[k][i] = t2
		M.Bfactors[k][j] = t1
	}
}

//Less returns true if the value in the Bfactors for
//atom i are less than that for atom j, and false otherwise.
func (M *Molecule) Less(i, j int) bool {
	return M.Bfactors[0][i] < M.Bfactors[0][j]
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
	if M != nil || M.Coords != nil || M.current < len(M.Coords) {

		return true
	}
	return false
}

//Next returns the  next frame and an error
func (M *Molecule) Next(V *v3.Matrix) error {
	if M.current >= len(M.Coords) {
		return newlastFrameError("", len(M.Coords)-1)
	}
	//	fmt.Println("CURR", M.current, len(M.Coords), V.NVecs(), M.Coords[M.current].NVecs()) ////////////////
	M.current++
	if V == nil {
		return nil
	}
	V.Copy(M.Coords[M.current-1])
	return nil
}

//InitRead initializes molecule to be read as a traj (not tested!)
func (M *Molecule) InitRead() error {
	if M == nil || len(M.Coords) == 0 {
		return CError{"Bad molecule", []string{"InitRead"}}
	}
	M.current = 0
	return nil
}

//NextConc takes a slice of bools and reads as many frames as elements the list has
//form the trajectory. The frames are discarted if the corresponding elemetn of the slice
//is false. The function returns a slice of channels through each of each of which
// a *matrix.DenseMatrix will be transmited
func (M *Molecule) NextConc(frames []*v3.Matrix) ([]chan *v3.Matrix, error) {
	toreturn := make([]chan *v3.Matrix, 0, len(frames))
	used := false

	for _, val := range frames {
		if val == nil {
			M.current++
			toreturn = append(toreturn, nil)
			continue
		}
		if M.current >= len(M.Coords) {
			lastframe := newlastFrameError("", len(M.Coords)-1)
			if used == false {
				return nil, lastframe
			} else {
				return toreturn, lastframe
			}
		}
		used = true
		toreturn = append(toreturn, make(chan *v3.Matrix))
		go func(a *v3.Matrix, pipe chan *v3.Matrix) {
			pipe <- a
		}(M.Coords[M.current], toreturn[len(toreturn)-1])
		M.current++
	}
	return toreturn, nil
}

//Close just sets the "curren" counter to 0.
//If you are using it as a trajectory, you can always just discard the molecule
//and let the CG take care of it, as there is nothing on disk linked to it..
func (M *Molecule) Close() {
	M.current = 0
}

/**End Traj interface implementation***********/

//End Molecule methods

//Traj Error

type lastFrameError struct {
	fileName string
	frame    int
	deco     []string
}

//Error returns an error message string.
func (E *lastFrameError) Error() string {
	return "EOF" //: Last frame in mol-based trajectory from file %10s reached at frame %10d", E.fileName, E.frame)
}

//Format returns the format used by the trajectory that returned the error.
func (E *lastFrameError) Format() string {
	return "mol"
}

//Frame returns the frame at which the error was detected.
func (E *lastFrameError) Frame() int {
	return E.frame
}

func (E *lastFrameError) Critical() bool {
	return false
}

//FileName returns the name of the file from where the trajectory that gave the error is read.
func (E *lastFrameError) FileName() string {
	return E.fileName
}

//NormalLastFrameTermination does nothing, it is there so we can have an interface unifying all
//"normal termination" errors so they can be filtered out by type switch.
func (E *lastFrameError) NormalLastFrameTermination() {
}

//Decorate will add the dec string to the decoration slice of strings of the error,
//and return the resulting slice.
func (E *lastFrameError) Decorate(dec string) []string {
	if dec == "" {
		return E.deco
	}
	E.deco = append(E.deco, dec)
	return E.deco
}

func newlastFrameError(filename string, frame int) *lastFrameError {
	e := new(lastFrameError)
	e.fileName = filename
	e.frame = frame
	return e
}

//End Traj Error

//The general concrete error type for the package

//CError (Concrete Error) is the concrete error type
//for the chem package, that implements chem.Error
type CError struct {
	msg  string
	deco []string
}

func (err CError) Error() string { return err.msg }

//Decorate will add the dec string to the decoration slice of strings of the error,
//and return the resulting slice.
func (err CError) Decorate(dec string) []string {
	if dec == "" {
		return err.deco
	}
	err.deco = append(err.deco, dec)
	return err.deco
}

//errDecorate will decorate a chem.Error, or use the message of the error
//and the name of the called to produce a new chem.Error (if the original error does not
//implement chem.Error
func errDecorate(err error, caller string) Error {
	if err == nil {
		return nil
	}
	err2, ok := err.(Error) //I know that is the type returned byt initRead
	if ok {
		err2.Decorate(caller)
		return err2
	}
	err3 := CError{err.Error(), []string{caller}}
	return err3
}

//PanicMsg is the type used for all the panics raised in the chem package
//so they can be easily recovered if needed. goChem only raises panics
//for programming errors: Attempts to to out of a matrice's bounds,
//dimension mismatches in matrices, etc.
type PanicMsg string

//Error returns a string with an error message
func (v PanicMsg) Error() string { return string(v) }

const (
	ErrNilData          = PanicMsg("goChem: Nil data given ")
	ErrInconsistentData = PanicMsg("goChem: Inconsistent data length ")
	ErrNilMatrix        = PanicMsg("goChem: Attempted to access nil v3.Matrix or gonum/mat64.Dense")
	ErrNilAtoms         = PanicMsg("goChem: Topology has a nil []*Atom slice")
	ErrNilAtom          = PanicMsg("goChem: Attempted to copy from or to a nil Atom")
	ErrAtomOutOfRange   = PanicMsg("goChem: Requested/Attempted setting Atom out of range")
	ErrNilFrame         = PanicMsg("goChem: Attempted to acces nil frame")
	ErrNotXx3Matrix     = PanicMsg("goChem: A v3.Matrix should have 3 columns")
	ErrCliffordRotation = PanicMsg("goChem-Clifford: Target and Result matrices must have the same dimensions. They cannot reference the same matrix") //the only panic that the Clifford functions throw.
)
