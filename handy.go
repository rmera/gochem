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
import "reflect"

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
		if isInInt(residues, at.Molid) && isInString(chains, string(at.Chain)) {
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

//IsIn returns the position of test in the slice set, or
// -1 if test is not present in set. Panics if set is not a slice
func isIn(test interface{}, set interface{}) int {
	vset := reflect.ValueOf(set)
	if reflect.TypeOf(set).Kind().String() != "slice" {
		panic("IsIn function needs a slice as second argument!")
	}
	if vset.Len() < 0 {
		return 1
	}
	for i := 0; i < vset.Len(); i++ {
		vcomp := vset.Index(i)
		comp := vcomp.Interface()
		if reflect.DeepEqual(test, comp) {
			return i
		}
	}
	return -1
}
