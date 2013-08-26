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

func Deg2Rad(f float64) float64 {
	return f * 0.0174533
}

func Rad2Deg(f float64) float64 {
	return f / 0.0174533
}

//Molecules2Atoms gets a selection list from a list of residues.
//It select all the atoms that form part of the residues in the list.
//It doesnt return errors, if a residue is out of range, no atom will
//be returned for it. Atoms are also required to be part of one of the chains
//specified in chains.
func Molecules2Atoms(mol Ref, residues []int, chains []string) []int {
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
func OnesMass(lenght int) *CoordMatrix {
	return gnOnes(lenght, 1)
}

//Super determines the best rotation and translations to superimpose the coords in test
//listed in testlst on te atoms of molecule templa, frame frametempla, listed in templalst.
//It applies those rotation and translations to the whole frame frametest of molecule test, in palce.
//testlst and templalst must have the same number of elements.
func Super(test, templa *CoordMatrix, testlst, templalst []int) (*CoordMatrix, error) {
	_, testcols := test.Dims()
	_, templacols := templa.Dims()

	ctest := gnZeros(len(testlst), testcols)
	if len(testlst) != 0 {
		ctest.SomeRows(test, testlst)
	}
	ctempla := gnZeros(len(templalst), templacols)
	if len(templalst) != 0 {
		ctempla.SomeRows(templa, templalst)
	}
	if len(templalst) != len(testlst) {
		return nil, fmt.Errorf("Mismatched template and test atom numbers: %d, %d", len(templalst), len(testlst))
	}
	_, rotation, trans1, trans2, err1 := GetSuper(ctest, ctempla)
	if err1 != nil {
		return nil, err1
	}
	test.AddRow(test, trans1)
	test = gnMul(test, rotation)
	test.AddRow(test, trans2)
	return test, nil
}

//Rotate about rotates the coordinates in coordsorig around by angle radians around the axis
//given by the vector axis. It returns the rotated coordsorig, since the original is not affected.
//Uses Clifford algebra.
func RotateAbout(coordsorig, ax1, ax2 *CoordMatrix, angle float64) (*CoordMatrix, error) {
	coords := gnClone(coordsorig)
	translation := gnClone(ax1)
	axis := gnClone(ax2)
	axis.Sub(axis, ax1) //now it became the rotation axis
	f := func() { coords.SubRow(coords, translation) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, err
	}
	Rot := Rotate(coords, axis, angle)
	g := func() { Rot.AddRow(Rot, translation) }
	if err := gnMaybe(gnPanicker(g)); err != nil {
		return nil, err
	}
	return Rot, nil
}

//EulerRotateAbout uses Euler angles to rotate the coordinates in coordsorig around by angle
//radians around the axis given by the vector axis. It returns the rotated coordsorig,
//since the original is not affected. It seems more clunky than the RotateAbout, which uses Clifford algebra.
//I leave it for benchmark, mostly, and might remove it later.
func EulerRotateAbout(coordsorig, ax1, ax2 *CoordMatrix, angle float64) (*CoordMatrix, error) {
	coords := gnClone(coordsorig)
	translation := gnClone(ax1)
	axis := gnClone(ax2)
	axis.Sub(axis, ax1) //now it became the rotation axis
	f := func() { coords.SubRow(coords, translation) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, err
	}
	Zswitch := GetSwitchZ(axis)
	coords = gnMul(coords, Zswitch) //rotated
	Zrot, err := GetRotateAroundZ(angle)
	if err != nil {
		return nil, err
	}
	Zsr, Zsc := Zswitch.Dims()
	RevZ := gnZeros(Zsr, Zsc)
	g := func() { RevZ.Inv(Zswitch) }
	if err := gnMaybe(gnPanicker(g)); err != nil {
		return nil, err
	}
	coords = gnMul(coords, Zrot) //rotated
	coords = gnMul(coords, RevZ)
	coords.AddRow(coords, translation)
	return coords, err
}

//Corrupted is a convenience function to check that a reference and a trajectory have the same number of atoms
func Corrupted(R Ref, X Traj) error {
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
func MakeWater(a1, a2 *CoordMatrix, distance, angle float64, oxygen bool) *CoordMatrix {
	water:=Zeros(3,3)
	WaterOHDist:=0.96
	WaterAngle:=52.25
	deg2rad := 0.0174533
	w:=EmptyCoords()
	w.RowView(water,0) //we first set the O coordinates
	w.Clone(a2)
	w.Sub(w,a1)
	w.Unit(w)
	dist:=Zeros(1,3)
	dist.Sub(a1,a2)
	a1a2dist:=dist.Norm(2)
	w.Scale(distance+a1a2dist,w)
	w.Add(w,a1)
	for i:=0;i<=1;i++{
		o:=EmptyCoords()
		o.RowView(water,0)
		w.RowView(water,i+1)
		w.Clone(o)
		w.Sub(w,a2)
		w.Unit(w)
		w.Scale(WaterOHDist+distance,w)
		o.Sub(o,a2)
		upp,_:=Cross3D(w,NewCoords([]float64{0,0,1},1,3))
		upp.Add(upp,o)
		upp.Add(upp,a2)
		//water.SetMatrix(3,0,upp)
		w.Add(w,a2)
		o.Add(o,a2)
		sign:=1.0
		if i==1{
			sign=-1.0
		}
		temp,_:=RotateAbout(w,o,upp,deg2rad*WaterAngle*sign)
		w.SetMatrix(0,0,temp)
	}
	var v1, v2 *CoordMatrix
	if angle!=0{
		v1=Zeros(1,3)
		v2=Zeros(1,3)
		v1.Sub(a2,a1)
		v2.Clone(v1)
		v2.Set(0,2,v2.At(0,2)+1) //a "random" modification. The idea is that its not colinear with v1
		v3,_:=Cross3D(v1,v2)
		v3.Add(v3,a2)
		water,_=RotateAbout(water,a2,v3,angle)
	}
	if oxygen{
		return water
		}
	//we move things so an hydrogen points to a2 and modify the distance acordingly.
	e1:=EmptyCoords()
	e2:=EmptyCoords()
	e3:=EmptyCoords()
	e1.RowView(water,0)
	e2.RowView(water,1)
	e3.RowView(water,2)
	if v1==nil{
		v1=Zeros(1,3)
		}
	if v2==nil{
		v2=Zeros(1,3)
		}
	v1.Sub(e2,e1)
	v2.Sub(e3,e1)
	axis,_:=Cross3D(v1,v2)
	axis.Add(axis,e1)
	water,_=RotateAbout(water,e1,axis,deg2rad*(180-WaterAngle))
	v1.Sub(e1,a2)
	v1.Unit(v1)
	v1.Scale(WaterOHDist,v1)
	water.AddRow(water,v1)
	return water
}


//This function will put the internal numbering+1 in the atoms and residue fields, so they match the current residues/atoms
//in the molecule
func FixNumbering(r Ref){
	resid:=0
	prevres:=-1
	for i:=0;i<r.Len();i++{
		at:=r.Atom(i)
		at.Id=i+1
		if prevres!=at.Molid{
			prevres=at.Molid
			resid++
		}
		at.Molid=resid
	}
}


//Takes a list of lists of residues and deletes from r
//all atoms not in the list or not belonging to the chain chain. 
//It caps the N and C terminal
//of each list with -COH for the N terminal and NH2 for C terminal.
//the residues on each sublist should be contiguous to each other.
//for instance, {6,7,8} is a valid sublist, {6,8,9} is not.
//This is NOT currently checked by the function!
func CutBackRef(r Ref, chain string, list [][]int) error{
	//i:=r.Len()
	for k,v:=range(list){
		nter:=v[0]
		cter:=v[len(v)-1]
		nresname:=""
		cresname:=""
		for j:=0;j<r.Len();j++{
			if r.Atom(j).Molid==nter && r.Atom(j).Chain==chain{
				nresname=r.Atom(j).Molname
				break
				}
		if nresname==""{
			//we will protest if the Nter is not found. If Cter is not found we will just
			//cut at the real Cter 
			return fmt.Errorf("list %d contains residue numbers out of boundaries",k)
			}
		}
		for j:=0;j<r.Len();j++{
			curr:=r.Atom(j)
			if curr.Chain!=chain{
				continue
				}
			if curr.Molid==cter{
				cresname=curr.Molname
				}
			if curr.Molid==nter-1{
				makeNcap(curr,nresname)
				}
			if curr.Molid==cter+1{
				makeCcap(curr,cresname)
				}
		}
	
		for j:=0;j<r.Len();j++{
			remove:=true
			for _,i:=range(list){
				if isInInt(i, r.Atom(j).Molid) && r.Atom(j).Chain==chain{
					remove=false
					break
					}
				}
			if remove{
				r.DelAtom(j)
				}
		}
	
	}
	return nil
}

func makeNcap(at *Atom,resname string){
	if  at.Name=="C" || at.Name=="O"{
		at.Molid=at.Molid+1
		at.Molname=resname
		} 
	if at.Name=="CA"{
		at.Name="HC"
		at.Symbol="H"
		at.Molid=at.Molid+1
		at.Molname=resname
		}
}

func makeCcap(at *Atom, resname string){
	if at.Name=="H" || at.Name=="N"{
		at.Molid=at.Molid-1
		at.Molname=resname
		}
	if at.Name=="CA"{
		at.Name="HN2"
		at.Symbol="H"
		at.Molid=at.Molid-1
		at.Molname=resname
		}
}


//Takes a list of lists of residues and produces a new set of coordinates
//whitout any atom not in the lists or not from the chain chain. It caps the N and C terminal
//of each list with -COH for the N terminal and NH2 for C terminal.
//the residues on each sublist should be contiguous to each other.
//for instance, {6,7,8} is a valid sublist, {6,8,9} is not.
//This is NOT currently checked by the function!
//In addition, the Ref provided should have already been processed by 
//CutBackRef, which is not checked either.
func CutBackCoords(r Ref, coords *CoordMatrix, chain string, list [][]int) (*CoordMatrix, error){
	//this is actually a really silly function. So far I dont check for errors, but I keep the return balue
	//In case I do later.
	var biglist []int
	for _,i:=range(list){
		smallist:=Molecules2Atoms(r,i,[]string{chain})
		biglist=append(biglist,smallist...)
	}
	newcoords:=Zeros(len(biglist),3)
	newcoords.SomeRows(coords,biglist)
	return newcoords, nil


}


