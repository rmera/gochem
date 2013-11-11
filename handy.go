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

//Molecules2Atoms gets a selection list from a list of residues.
//It select all the atoms that form part of the residues in the list.
//It doesnt return errors, if a residue is out of range, no atom will
//be returned for it. Atoms are also required to be part of one of the chains
//specified in chains.
func Molecules2Atoms(mol Atomer, residues []int, chains []string) []int {
	atlist := make([]int, 0, len(residues)*3)
	for key := 0; key < mol.Len(); key++ {
		at := mol.Atom(key)
		if isInInt(residues, at.Molid) && (isInString(chains, at.Chain) || len(chains) == 0) {
			atlist = append(atlist, key)
		}
	}
	return atlist

}

//Ones mass returns a column matrix with lenght rosw.
//This matrix can be used as a dummy mass matrix
//for geometric calculations.
func OnesMass(lenght int) *VecMatrix {
	return chemDense2VecMatrix(gnOnes(lenght, 1))
}

//Super determines the best rotation and translations to superimpose the coords in test
//listed in testlst on te atoms of molecule templa, frame frametempla, listed in templalst.
//It applies those rotation and translations to the whole frame frametest of molecule test, in palce.
//testlst and templalst must have the same number of elements.
func Super(test, templa *VecMatrix, testlst, templalst []int) (*VecMatrix, error) {
	//_, testcols := test.Dims()
	//_, templacols := templa.Dims()
	if len(testlst) != len(templalst) {
		return nil, fmt.Errorf("testlst and templalst must have the same lenght. testlst: %d, templalst %d", len(testlst), len(templalst))
	}
	ctest := ZeroVecs(len(testlst))
	if len(testlst) != 0 {
		ctest.SomeVecs(test, testlst)
	}
	ctempla := ZeroVecs(len(templalst))
	if len(templalst) != 0 {
		ctempla.SomeVecs(templa, templalst)
	}
	if len(templalst) != len(testlst) {
		return nil, fmt.Errorf("Mismatched template and test atom numbers: %d, %d", len(templalst), len(testlst))
	}
	_, rotation, trans1, trans2, err1 := RotatorTranslatorToSuper(ctest, ctempla)
	if err1 != nil {
		return nil, err1
	}
	test.AddVec(test, trans1)
	//	fmt.Println("test1",test, rotation) /////////////77
	test.Mul(test, rotation)
	//	fmt.Println("test2",test) ///////////
	test.AddVec(test, trans2)
	//	fmt.Println("test3",test) ///////
	return test, nil
}

//Rotate about rotates the coordinates in coordsorig around by angle radians around the axis
//given by the vector axis. It returns the rotated coordsorig, since the original is not affected.
//Uses Clifford algebra.
func RotateAbout(coordsorig, ax1, ax2 *VecMatrix, angle float64) (*VecMatrix, error) {
	coords := ZeroVecs(coordsorig.NVecs())
	translation := ZeroVecs(ax1.NVecs())
	translation.Copy(ax1)
	axis := ZeroVecs(ax2.NVecs())
	axis.Sub(ax2, ax1) // the rotation axis
	f := func() { coords.SubVec(coordsorig, translation) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, err
	}
	Rot := Rotate(coords, axis, angle)
	g := func() { Rot.AddVec(Rot, translation) }
	if err := gnMaybe(gnPanicker(g)); err != nil {
		return nil, err
	}
	return Rot, nil
}

//EulerRotateAbout uses Euler angles to rotate the coordinates in coordsorig around by angle
//radians around the axis given by the vector axis. It returns the rotated coordsorig,
//since the original is not affected. It seems more clunky than the RotateAbout, which uses Clifford algebra.
//I leave it for benchmark, mostly, and might remove it later. There is no test for this function!
func EulerRotateAbout(coordsorig, ax1, ax2 *VecMatrix, angle float64) (*VecMatrix, error) {
	r, _ := coordsorig.Dims()
	coords := ZeroVecs(r)
	translation := ZeroVecs(ax1.NVecs())
	translation.Copy(ax1)
	axis := ZeroVecs(ax2.NVecs())
	axis.Sub(ax2, ax1) //now it became the rotation axis
	f := func() { coords.SubVec(coordsorig, translation) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, err
	}
	Zswitch := RotatorToNewZ(axis)
	coords.Mul(coords, Zswitch) //rotated
	Zrot, err := RotatorAroundZ(angle)
	if err != nil {
		return nil, err
	}
	Zsr, _ := Zswitch.Dims()
	RevZ := ZeroVecs(Zsr)
	g := func() { RevZ = gnInverse(Zswitch) }
	if err := gnMaybe(gnPanicker(g)); err != nil {
		return nil, err
	}
	coords.Mul(coords, Zrot) //rotated
	coords.Mul(coords, RevZ)
	coords.AddVec(coords, translation)
	return coords, err
}

//Corrupted is a convenience function to check that a reference and a trajectory have the same number of atoms
func Corrupted(X Traj, R Atomer) error {
	if X.Len() != R.Len() {
		return fmt.Errorf("Mismatched number of atoms/coordinates")
	}
	return nil
}

//Some internal convenience functions.

//isIn is a helper for the RamaList function,
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

//Creates a water molecule at distance Angstroms from a2, in a direction that is angle radians from the axis defined by a1 and a2.
//Notice that the exact position of the water is not well defined when angle is not zero. One can always use the RotateAbout
//function to move the molecule to the desired location. If oxygen is true, the oxygen will be pointing to a2. Otherwise,
//one of the hydrogens will.
func MakeWater(a1, a2 *VecMatrix, distance, angle float64, oxygen bool) *VecMatrix {
	water := ZeroVecs(3)
	WaterOHDist := 0.96
	WaterAngle := 52.25
	deg2rad := 0.0174533
	w := water.VecView(0) //we first set the O coordinates
	w.Copy(a2)
	w.Sub(w, a1)
	w.Unit(w)
	dist := ZeroVecs(1)
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
		t, _ := NewVecs([]float64{0, 0, 1})
		upp := ZeroVecs(1)
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
	var v1, v2 *VecMatrix
	if angle != 0 {
		v1 = ZeroVecs(1)
		v2 = ZeroVecs(1)
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
		v1 = ZeroVecs(1)
	}
	if v2 == nil {
		v2 = ZeroVecs(1)
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

//This function will put the internal numbering+1 in the atoms and residue fields, so they match the current residues/atoms
//in the molecule
func FixNumbering(r Ref) {
	resid := 0
	prevres := -1
	for i := 0; i < r.Len(); i++ {
		at := r.Atom(i)
		at.Id = i + 1
		if prevres != at.Molid {
			prevres = at.Molid
			resid++
		}
		at.Molid = resid
	}
}

//Takes a list of lists of residues and deletes from r
//all atoms not in the list or not belonging to the chain chain.
//It caps the N and C terminal
//of each list with -COH for the N terminal and NH2 for C terminal.
//the residues on each sublist should be contiguous to each other.
//for instance, {6,7,8} is a valid sublist, {6,8,9} is not.
//This is NOT currently checked by the function!. It returns the list of kept atoms
func CutBackRef(r Atomer, chains []string, list [][]int) ([]int, error) {
	//i:=r.Len()
	if len(chains) != len(list) {
		return nil, fmt.Errorf("Mismatched chains (%d) and list (%d) slices", len(chains), len(list))
	}
	var ret []int //This will be filled with the atoms that are kept, and will be returned.
	for k, v := range list {
		nter := v[0]
		cter := v[len(v)-1]
		nresname := ""
		cresname := ""
		for j := 0; j < r.Len(); j++ {
			if r.Atom(j).Molid == nter && r.Atom(j).Chain == chains[k] {
				nresname = r.Atom(j).Molname
				break
			}
		}
		if nresname == "" {
			//we will protest if the Nter is not found. If Cter is not found we will just
			//cut at the real Cter
			return nil, fmt.Errorf("list %d contains residue numbers out of boundaries", k)
		}
		for j := 0; j < r.Len(); j++ {
			curr := r.Atom(j)
			if curr.Chain != chains[k] {
				continue
			}
			if curr.Molid == cter {
				cresname = curr.Molname
			}
			if curr.Molid == nter-1 {
				makeNcap(curr, nresname)
			}
			if curr.Molid == cter+1 {
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
	at.Molid = at.Molid + 1
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
	at.Molid = at.Molid - 1
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
func CutBackCoords(r Ref, coords *VecMatrix, chain string, list [][]int) (*VecMatrix, error) {
	//this is actually a really silly function. So far I dont check for errors, but I keep the return balue
	//In case I do later.
	var biglist []int
	for _, i := range list {
		smallist := Molecules2Atoms(r, i, []string{chain})
		biglist = append(biglist, smallist...)
	}
	NewVecs := ZeroVecs(len(biglist), 3)
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
		if isInInt(list, curr.Molid) && isInString(chain, curr.Chain) {
			if curr.Name == "CB" {
				//			pairs[len(pairs)-1][1]=i //I am assuming that CA will show before CB in the PDB, which is rather weak
				//		paairs=append(pairs,make([]int,1,2))
			}
			if curr.Name == "CA" {
				curr.Name = "HB4"
				curr.Symbol = "H"
				//		pairs[len(pairs)-1]=append(pairs[len(pairs)-1],i)
			} else if isInString([]string{"C", "H", "HA", "O", "N"}, curr.Name) { //change the res number of the backbone so it is not considered
				curr.Molid = -1
			}

		}
	}
	newlist := Molecules2Atoms(r, list, chain)
	return newlist
}

func CutAlphaRef(r Atomer, chain []string, list []int) []int {
	for i := 0; i < r.Len(); i++ {
		curr := r.Atom(i)
		if isInInt(list, curr.Molid) && isInString(chain, curr.Chain) {
			if curr.Name == "C" {
				curr.Name = "HA2"
				curr.Symbol = "H"
			} else if curr.Name == "N" {
				curr.Name = "HA3"
				curr.Symbol = "H"
			} else if isInString([]string{"H", "O"}, curr.Name) { //change the res number of the backbone so it is not considered
				curr.Molid = -1
			}

		}
	}
	newlist := Molecules2Atoms(r, list, chain)
	return newlist
}

//This will tag all atoms with a given name in a given list of atoms.
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

//Scales all bonds between atoms in the same residue with names n1, n2 to a final lenght finallengt, by moving the atoms n2.
//the operation is executed in place.
func ScaleBonds(coords *VecMatrix, mol Atomer, n1, n2 string, finallenght float64) {
	for i := 0; i < mol.Len(); i++ {
		c1 := mol.Atom(i)
		if c1.Name != n1 {
			continue
		}
		for j := 0; j < mol.Len(); j++ {
			c2 := mol.Atom(j)
			if c1.Molid == c2.Molid && c1.Name == n1 && c2.Name == n2 {
				A := coords.VecView(i)
				B := coords.VecView(j)
				ScaleBond(A, B, finallenght)
			}
		}
	}
}

//ScaleCHBond takes a bond and moves the H (in place) so the distance between them is bond.
//CAUTION: I have only tested it for the case where the original distance>bond, although I think it will also work in the other case.
func ScaleBond(C, H *VecMatrix, bond float64) {
	Odist := ZeroVecs(1)
	Odist.Sub(H, C)
	distance := Odist.Norm(0)
	Odist.Scale((distance-bond)/distance, Odist)
	H.Sub(H, Odist)
}

func MergeAtomers(A, B Atomer) (*Topology, error) {
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
	a, aok := A.(ReadRef)
	b, bok := B.(ReadRef)
	var charge, multi int
	if aok && bok {
		charge = a.Charge() + b.Charge()
		multi = (a.Multi() - 1 + b.Multi()) //Not TOO sure about this.
	} else {
		multi = 1
	}
	r, err := NewTopology(full, charge, multi)
	return r, err
}
