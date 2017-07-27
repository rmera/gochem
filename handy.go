/*
 * files.go, part of gochem.
 *
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
 *
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package chem

import "fmt"
import "math"
import "github.com/rmera/gochem/v3"

//Molecules2Atoms gets a selection list from a list of residues.
//It select all the atoms that form part of the residues in the list.
//It doesnt return errors, if a residue is out of range, no atom will
//be returned for it. Atoms are also required to be part of one of the chains
//specified in chains.
func Molecules2Atoms(mol Atomer, residues []int, chains []string) []int {
	atlist := make([]int, 0, len(residues)*3)
	for key := 0; key < mol.Len(); key++ {
		at := mol.Atom(key)
		if isInInt(residues, at.MolID) && (isInString(chains, at.Chain) || len(chains) == 0) {
			atlist = append(atlist, key)
		}
	}
	return atlist

}

//MolIDNameChain2Index takes a molID (residue number), atom name, chain index and a molecule Atomer.
//it returns the index associated with the atom in question in the Ref. The function returns also an error (if failure of warning)
// or nil (if succses and no warnings). Note that this function is not efficient to call several times to retrieve many atoms.
func MolIDNameChain2Index(mol Atomer, molID int, name, chain string) (int, error) {
	var ret int = -1
	var err error
	if mol == nil {
		return -1, CError{"goChem: Given a nil chem.Atomer", []string{"MolIDNameChain2Index"}}
	}
	for i := 0; i != mol.Len(); i++ {
		a := mol.Atom(i)
		if a.Name == "" && err == nil {
			err = CError{"Warning: The Atoms does not seem to contain PDB-type information", []string{"MolIDNameChain2Index"}} //We set this error but will still keep running the function in case the data is present later in the molecule.
		}
		if a.MolID == molID && a.Name == name && a.Chain == chain {
			ret = i
			break
		}

	}
	if ret == -1 {
		var p string
		if err != nil {
			p = err.Error()
		}
		err = CError{fmt.Sprintf("%s.  No atomic index found in the Atomer given for the given MolID, atom name and chain.", p), []string{"MolIDNameChain2Index"}}
	}
	return ret, err
}

//OnesMass returns a column matrix with lenght rosw.
//This matrix can be used as a dummy mass matrix
//for geometric calculations.
func OnesMass(lenght int) *v3.Matrix {
	return v3.Dense2Matrix(gnOnes(lenght, 1))
}

//Super determines the best rotation and translations to superimpose the coords in test
//considering only the atoms present in the slices of int slices indexes.
//The first indexes slices will be assumed to contain test indexes and the second, template indexes.
//If you give only one, it will be assumed to correspondo to whatever molecule
//that has more atoms than the elements in the slice. The same number of atoms
//has to be considered for superposition in both systems.
//It applies those rotation and translations to the whole molecule test, in palce.
//testlst and templalst must have the same number of elements. 
func Super(test, templa *v3.Matrix, indexes ...[]int) (*v3.Matrix, error) {
	var ctest *v3.Matrix
	var ctempla *v3.Matrix
	if len(indexes)==0{
		ctest=test
		ctempla=templa
	}else if len(indexes)==1{
		if  test.NVecs()>len(indexes[0]){
			ctest = v3.Zeros(len(indexes[0]))
			ctest.SomeVecs(test, indexes[0])
			ctempla=templa
		}else if templa.NVecs()>len(indexes[0]){
			ctempla = v3.Zeros(len(indexes[0]))
			ctempla.SomeVecs(templa, indexes[0])
		}else{
			return nil, fmt.Errorf("chem.Super: Indexes don't match molecules")
		}
	}else{
		ctest = v3.Zeros(len(indexes[0]))
		ctest.SomeVecs(test, indexes[0])
		ctempla = v3.Zeros(len(indexes[1]))
		ctempla.SomeVecs(templa, indexes[1])
	}

	if ctest.NVecs() != ctempla.NVecs() {
		return nil, fmt.Errorf("chem.Super: Ill formed coordinates for Superposition")
	}

	_, rotation, trans1, trans2, err1 := RotatorTranslatorToSuper(ctest, ctempla)
	if err1 != nil {
		return nil, errDecorate(err1, "Super")
	}
	test.AddVec(test, trans1)
	//	fmt.Println("test1",test, rotation) /////////////77
	test.Mul(test, rotation)
	//	fmt.Println("test2",test) ///////////
	test.AddVec(test, trans2)
	//	fmt.Println("test3",test) ///////
	return test, nil
}

//RotateAbout about rotates the coordinates in coordsorig around by angle radians around the axis
//given by the vector axis. It returns the rotated coordsorig, since the original is not affected.
//Uses Clifford algebra.
func RotateAbout(coordsorig, ax1, ax2 *v3.Matrix, angle float64) (*v3.Matrix, error) {
	coordsLen := coordsorig.NVecs()
	coords := v3.Zeros(coordsLen)
	translation := v3.Zeros(ax1.NVecs())
	translation.Copy(ax1)
	axis := v3.Zeros(ax2.NVecs())
	axis.Sub(ax2, ax1) // the rotation axis
	f := func() { coords.SubVec(coordsorig, translation) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, CError{err.Error(), []string{"v3.Matrix.SubVec", "RotateAbout"}}
	}
	Rot := v3.Zeros(coordsLen)
	Rot = Rotate(coords, Rot, axis, angle)
	g := func() { Rot.AddVec(Rot, translation) }
	if err := gnMaybe(gnPanicker(g)); err != nil {
		return nil, CError{err.Error(), []string{"v3.Matrix.AddVec", "RotateAbout"}}

	}
	return Rot, nil
}

//EulerRotateAbout uses Euler angles to rotate the coordinates in coordsorig around by angle
//radians around the axis given by the vector axis. It returns the rotated coordsorig,
//since the original is not affected. It seems more clunky than the RotateAbout, which uses Clifford algebra.
//I leave it for benchmark, mostly, and might remove it later. There is no test for this function!
func EulerRotateAbout(coordsorig, ax1, ax2 *v3.Matrix, angle float64) (*v3.Matrix, error) {
	r, _ := coordsorig.Dims()
	coords := v3.Zeros(r)
	translation := v3.Zeros(ax1.NVecs())
	translation.Copy(ax1)
	axis := v3.Zeros(ax2.NVecs())
	axis.Sub(ax2, ax1) //now it became the rotation axis
	f := func() { coords.SubVec(coordsorig, translation) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, CError{err.Error(), []string{"v3.Matrix.Subvec", "EulerRotateAbout"}}

	}
	Zswitch := RotatorToNewZ(axis)
	coords.Mul(coords, Zswitch) //rotated
	Zrot, err := RotatorAroundZ(angle)
	if err != nil {
		return nil, errDecorate(err, "EulerRotateAbout")
	}
	//	Zsr, _ := Zswitch.Dims()
	//	RevZ := v3.Zeros(Zsr)
	RevZ, err := gnInverse(Zswitch, nil)
	if err != nil {
		return nil, errDecorate(err, "EulerRotateAbout")
	}
	coords.Mul(coords, Zrot) //rotated
	coords.Mul(coords, RevZ)
	coords.AddVec(coords, translation)
	return coords, nil
}

//Corrupted is a convenience function to check that a reference and a trajectory have the same number of atoms
func Corrupted(X Traj, R Atomer) error {
	if X.Len() != R.Len() {
		return CError{"Mismatched number of atoms/coordinates", []string{"Corrupted"}}
	}
	return nil
}

//Some internal convenience functions.

//isInInt is a helper for the RamaList function,
//returns true if test is in container, false otherwise.
func isInInt(container []int, test int) bool {
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
func isInString(container []string, test string) bool {
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

//MakeWater Creates a water molecule at distance Angstroms from a2, in a direction that is angle radians from the axis defined by a1 and a2.
//Notice that the exact position of the water is not well defined when angle is not zero. One can always use the RotateAbout
//function to move the molecule to the desired location. If oxygen is true, the oxygen will be pointing to a2. Otherwise,
//one of the hydrogens will.
func MakeWater(a1, a2 *v3.Matrix, distance, angle float64, oxygen bool) *v3.Matrix {
	water := v3.Zeros(3)
	const WaterOHDist = 0.96
	const WaterAngle = 52.25
	const deg2rad = 0.0174533
	w := water.VecView(0) //we first set the O coordinates
	w.Copy(a2)
	w.Sub(w, a1)
	w.Unit(w)
	dist := v3.Zeros(1)
	dist.Sub(a1, a2)
	a1a2dist := dist.Norm(0)
	fmt.Println("ala2dist", a1a2dist, distance) ////////////////7777
	w.Scale(distance+a1a2dist, w)
	w.Add(w, a1)
	for i := 0; i <= 1; i++ {
		o := water.VecView(0)
		w = water.VecView(i + 1)
		w.Copy(o)
		fmt.Println("w1", w) ////////
		w.Sub(w, a2)
		fmt.Println("w12", w) ///////////////
		w.Unit(w)
		fmt.Println("w4", w)
		w.Scale(WaterOHDist+distance, w)
		fmt.Println("w3", w, WaterOHDist, distance)
		o.Sub(o, a2)
		t, _ := v3.NewMatrix([]float64{0, 0, 1})
		upp := v3.Zeros(1)
		upp.Cross(w, t)
		fmt.Println("upp", upp, w, t)
		upp.Add(upp, o)
		upp.Add(upp, a2)
		//water.SetMatrix(3,0,upp)
		w.Add(w, a2)
		o.Add(o, a2)
		sign := 1.0
		if i == 1 {
			sign = -1.0
		}
		temp, _ := RotateAbout(w, o, upp, deg2rad*WaterAngle*sign)
		w.SetMatrix(0, 0, temp)
	}
	var v1, v2 *v3.Matrix
	if angle != 0 {
		v1 = v3.Zeros(1)
		v2 = v3.Zeros(1)
		v1.Sub(a2, a1)
		v2.Copy(v1)
		v2.Set(0, 2, v2.At(0, 2)+1) //a "random" modification. The idea is that its not colinear with v1
		v3 := cross(v1, v2)
		v3.Add(v3, a2)
		water, _ = RotateAbout(water, a2, v3, angle)
	}
	if oxygen {
		return water
	}
	//we move things so an hydrogen points to a2 and modify the distance acordingly.
	e1 := water.VecView(0)
	e2 := water.VecView(1)
	e3 := water.VecView(2)
	if v1 == nil {
		v1 = v3.Zeros(1)
	}
	if v2 == nil {
		v2 = v3.Zeros(1)
	}
	v1.Sub(e2, e1)
	v2.Sub(e3, e1)
	axis := cross(v1, v2)
	axis.Add(axis, e1)
	water, _ = RotateAbout(water, e1, axis, deg2rad*(180-WaterAngle))
	v1.Sub(e1, a2)
	v1.Unit(v1)
	v1.Scale(WaterOHDist, v1)
	water.AddVec(water, v1)
	return water
}

//FixNumbering will put the internal numbering+1 in the atoms and residue fields, so they match the current residues/atoms
//in the molecule
func FixNumbering(r Atomer) {
	resid := 0
	prevres := -1
	for i := 0; i < r.Len(); i++ {
		at := r.Atom(i)
		at.ID = i + 1
		if prevres != at.MolID {
			prevres = at.MolID
			resid++
		}
		at.MolID = resid
	}
}

//CutBackRef takes a list of lists of residues and deletes from r
//all atoms not in the list or not belonging to the chain chain.
//It caps the N and C terminal
//of each list with -COH for the N terminal and NH2 for C terminal.
//the residues on each sublist should be contiguous to each other.
//for instance, {6,7,8} is a valid sublist, {6,8,9} is not.
//This is NOT currently checked by the function!. It returns the list of kept atoms
func CutBackRef(r Atomer, chains []string, list [][]int) ([]int, error) {
	//i:=r.Len()
	if len(chains) != len(list) {
		return nil, CError{fmt.Sprintf("Mismatched chains (%d) and list (%d) slices", len(chains), len(list)), []string{"CutBackRef"}}
	}
	var ret []int //This will be filled with the atoms that are kept, and will be returned.
	for k, v := range list {
		nter := v[0]
		cter := v[len(v)-1]
		nresname := ""
		cresname := ""
		for j := 0; j < r.Len(); j++ {
			if r.Atom(j).MolID == nter && r.Atom(j).Chain == chains[k] {
				nresname = r.Atom(j).Molname
				break
			}
		}
		if nresname == "" {
			//we will protest if the Nter is not found. If Cter is not found we will just
			//cut at the real Cter
			return nil, CError{fmt.Sprintf("list %d contains residue numbers out of boundaries", k), []string{"CutBackRef"}}

		}
		for j := 0; j < r.Len(); j++ {
			curr := r.Atom(j)
			if curr.Chain != chains[k] {
				continue
			}
			if curr.MolID == cter {
				cresname = curr.Molname
			}
			if curr.MolID == nter-1 {
				makeNcap(curr, nresname)
			}
			if curr.MolID == cter+1 {
				makeCcap(curr, cresname)
			}
		}
	}
	for _, i := range list {
		t := Molecules2Atoms(r, i, chains)
		//	fmt.Println("t", len(t))
		ret = append(ret, t...)
	}
	//	j:=0
	//	for i:=0;;i++{
	//		index:=i-j
	//		if index>=r.Len(){
	//			break
	//		}
	//		if !isInInt(ret, i){
	//			r.DelAtom(index)
	//			j++
	//		}
	//	}
	return ret, nil
}

func makeNcap(at *Atom, resname string) {
	if !isInString([]string{"C", "O", "CA"}, at.Name) {
		return
	}
	at.MolID = at.MolID + 1
	at.Molname = resname
	if at.Name == "C" {
		at.Name = "CTZ"
	}
	if at.Name == "CA" {
		at.Name = "HCZ"
		at.Symbol = "H"
	}
}

func makeCcap(at *Atom, resname string) {
	if !isInString([]string{"N", "H", "CA"}, at.Name) {
		return
	}
	at.MolID = at.MolID - 1
	at.Molname = resname
	if at.Name == "N" {
		at.Name = "NTZ"
	}
	if at.Name == "CA" {
		at.Name = "HNZ"
		at.Symbol = "H"
	}
}

/*
//Takes a list of lists of residues and produces a new set of coordinates
//whitout any atom not in the lists or not from the chain chain. It caps the N and C terminal
//of each list with -COH for the N terminal and NH2 for C terminal.
//the residues on each sublist should be contiguous to each other.
//for instance, {6,7,8} is a valid sublist, {6,8,9} is not.
//This is NOT currently checked by the function!
//In addition, the Ref provided should have already been processed by
//CutBackRef, which is not checked either.
func CutBackCoords(r Ref, coords *v3.Matrix, chain string, list [][]int) (*v3.Matrix, error) {
	//this is actually a really silly function. So far I dont check for errors, but I keep the return balue
	//In case I do later.
	var biglist []int
	for _, i := range list {
		smallist := Molecules2Atoms(r, i, []string{chain})
		biglist = append(biglist, smallist...)
	}
	NewVecs := v3.Zeros(len(biglist), 3)
	NewVecs.SomeVecs(coords, biglist)
	return NewVecs, nil

}
*/

//CutLateralRef will return a list with the atom indexes of the lateral chains of the residues in list
//for each of these residues it will change the alpha carbon to oxygen and change the residue number of the rest
//of the backbone to -1.
func CutBetaRef(r Atomer, chain []string, list []int) []int {
	//	pairs := make([][]int,1,10)
	//	pairs[0]=make([]int,0,2)
	for i := 0; i < r.Len(); i++ {
		curr := r.Atom(i)
		if isInInt(list, curr.MolID) && isInString(chain, curr.Chain) {
			if curr.Name == "CB" {
				//			pairs[len(pairs)-1][1]=i //I am assuming that CA will show before CB in the PDB, which is rather weak
				//		paairs=append(pairs,make([]int,1,2))
			}
			if curr.Name == "CA" {
				curr.Name = "HB4"
				curr.Symbol = "H"
				//		pairs[len(pairs)-1]=append(pairs[len(pairs)-1],i)
			} else if isInString([]string{"C", "H", "HA", "O", "N"}, curr.Name) { //change the res number of the backbone so it is not considered
				curr.MolID = -1
			}

		}
	}
	newlist := Molecules2Atoms(r, list, chain)
	return newlist
}

//CutAlphaRef will return a list with the atoms in the residues indicated by list, in the chains given.
//The carbonyl carbon and amide nitrogen for each residue will be transformer into hydrogens. The MolID of the
//other backbone atoms will be set to -1 so they are no longer considered.
func CutAlphaRef(r Atomer, chain []string, list []int) []int {
	for i := 0; i < r.Len(); i++ {
		curr := r.Atom(i)
		if isInInt(list, curr.MolID) && isInString(chain, curr.Chain) {
			if curr.Name == "C" {
				curr.Name = "HA2"
				curr.Symbol = "H"
			} else if curr.Name == "N" {
				curr.Name = "HA3"
				curr.Symbol = "H"
			} else if isInString([]string{"H", "O"}, curr.Name) { //change the res number of the backbone so it is not considered
				curr.MolID = -1
			}

		}
	}
	newlist := Molecules2Atoms(r, list, chain)
	return newlist
}

//TagAtomsByName will tag all atoms with a given name in a given list of atoms.
//return the number of tagged atoms
func TagAtomsByName(r Atomer, name string, list []int) int {
	tag := 0
	for i := 0; i < r.Len(); i++ {
		curr := r.Atom(i)
		if isInInt(list, i) && curr.Name == name {
			curr.Tag = 1
			tag++
		}
	}
	return tag
}

//ScaleBonds scales all bonds between atoms in the same residue with names n1, n2 to a final lenght finallengt, by moving the atoms n2.
//the operation is executed in place.
func ScaleBonds(coords *v3.Matrix, mol Atomer, n1, n2 string, finallenght float64) {
	for i := 0; i < mol.Len(); i++ {
		c1 := mol.Atom(i)
		if c1.Name != n1 {
			continue
		}
		for j := 0; j < mol.Len(); j++ {
			c2 := mol.Atom(j)
			if c1.MolID == c2.MolID && c1.Name == n1 && c2.Name == n2 {
				A := coords.VecView(i)
				B := coords.VecView(j)
				ScaleBond(A, B, finallenght)
			}
		}
	}
}

//ScaleBond takes a C-H bond and moves the H (in place) so the distance between them is the one given (bond).
//CAUTION: I have only tested it for the case where the original distance>bond, although I expect it to also work in the other case.
func ScaleBond(C, H *v3.Matrix, bond float64) {
	Odist := v3.Zeros(1)
	Odist.Sub(H, C)
	distance := Odist.Norm(0)
	Odist.Scale((distance-bond)/distance, Odist)
	H.Sub(H, Odist)
}

//Merges A and B in a single topology which is returned
func MergeAtomers(A, B Atomer) *Topology {
	al := A.Len()
	l := al + B.Len()
	full := make([]*Atom, l, l)
	for k, _ := range full {
		if k < al {
			full[k] = A.Atom(k)
		} else {
			full[k] = B.Atom(k - al)
		}
	}
	a, aok := A.(AtomMultiCharger)
	b, bok := B.(AtomMultiCharger)
	var charge, multi int
	if aok && bok {
		charge = a.Charge() + b.Charge()
		multi = (a.Multi() - 1 + b.Multi()) //Not TOO sure about this.
	} else {
		multi = 1
	}
	return NewTopology(full, charge, multi)
}

//SelCone, Given a set of cartesian points in sellist, obtains a vector "plane" normal to the best plane passing through the points.
//It selects atoms from the set A that are inside a cone in the direction of "plane" that starts from the geometric center of the cartesian points,
//and has an angle of angle (radians), up to a distance distance. The cone is approximated by a set of radius-increasing cilinders with height thickness.
//If one starts from one given point, 2 cgnOnes, one in each direction, are possible. If whatcone is 0, both cgnOnes are considered.
//if whatcone<0, only the cone opposite to the plane vector direction. If whatcone>0, only the cone in the plane vector direction.
//the 'initial' argument  allows the construction of a truncate cone with a radius of initial.
func SelCone(B, selection *v3.Matrix, angle, distance, thickness, initial float64, whatcone int) []int {
	A := v3.Zeros(B.NVecs())
	A.Copy(B) //We will be altering the input so its better to work with a copy.
	ar, _ := A.Dims()
	selected := make([]int, 0, 3)
	neverselected := make([]int, 0, 30000)     //waters that are too far to ever be selected
	nevercutoff := distance / math.Cos(angle)  //cutoff to be added to neverselected
	A, _, err := MassCenter(A, selection, nil) //Centrate A in the geometric center of the selection, Its easier for the following calculations
	if err != nil {
		panic(PanicMsg(err.Error()))
	}
	selection, _, _ = MassCenter(selection, selection, nil) //Centrate the selection as well
	plane, err := BestPlane(selection, nil)                 //I have NO idea which direction will this vector point. We might need its negative.
	if err != nil {
		panic(PanicMsg(err.Error()))
	}
	for i := thickness / 2; i <= distance; i += thickness {
		maxdist := math.Tan(angle)*i + initial //this should give me the radius of the cone at this point
		for j := 0; j < ar; j++ {
			if isInInt(selected, j) || isInInt(neverselected, j) { //we dont scan things that we have already selected, or are too far
				continue
			}
			atom := A.VecView(j)
			proj := Projection(atom, plane)
			norm := proj.Norm(0)
			//Now at what side of the plane is the atom?
			angle := Angle(atom, plane)
			if whatcone > 0 {
				if angle > math.Pi/2 {
					continue
				}
			} else if whatcone < 0 {
				if angle < math.Pi/2 {
					continue
				}
			}
			if norm > i+(thickness/2.0) || norm < (i-thickness/2.0) {
				continue
			}
			proj.Sub(proj, atom)
			projnorm := proj.Norm(0)
			if projnorm <= maxdist {
				selected = append(selected, j)
			}
			if projnorm >= nevercutoff {
				neverselected = append(neverselected, j)
			}
		}
	}
	return selected
}
