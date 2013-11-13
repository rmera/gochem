/*
 * gochem_test.go
 *
 * Copyright 2013  <rmera@Holmes>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

package chem

import "fmt"
import "os"
import "testing"

//TestMultiXYZ tests that multi-XYZ files are opened and read correctly.
func TestXYZIO(Te *testing.T) {
	mol, err := XYZRead("test/sample.xyz")
	if err != nil {
		fmt.Println("There was an error!", err.Error())
		Te.Error(err)
	}
	fmt.Println("XYZ read!")
	XYZWrite("test/sampleFirst.xyz", mol.Coords[0], mol)
}

func TestPDBIO(Te *testing.T) {
	mol, err := PDBRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	err = PDBWrite("test/2c9vIO.pdb", mol.Coords[0], mol, mol.Bfactors[0])
	if err != nil {
		Te.Error(err)
	}
	//for the 3 residue  I should get -131.99, 152.49.
}

//TestChangeAxis reads the PDB 2c9v.pdb from the test directory, collects
//The CA and CB of residue D124 of the chain A, and uses Clifford algebra to rotate the
//whole molecule such as the vector defined by these 2 atoms is
//aligned with the Z axis. The new molecule is written
//as 2c9v_aligned.pdb to the test folder.
func TestChangeAxis(Te *testing.T) {
	mol, err := PDBRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	PDBWrite("test/2c9v-Readtest.pdb", mol.Coords[0], mol, nil)
	//The selection thing
	orient_atoms := [2]int{0, 0}
	for index := 0; index < mol.Len(); index++ {
		atom := mol.Atom(index)
		if atom.Chain == "A" && atom.Molid == 124 {
			if atom.Name == "CA" {
				orient_atoms[0] = index
			} else if atom.Name == "CB" {
				orient_atoms[1] = index
			}
		}
	}
	//Get the axis of rotation
	//ov1:=mol.Coord(orient_atoms[0], 0)
	ov2 := mol.Coord(orient_atoms[1], 0)
	//now we center the thing in the beta carbon of D124
	mol.Coords[0].SubVec(mol.Coords[0], ov2)
	PDBWrite("test/2c9v-translated.pdb", mol.Coords[0], mol, nil)
	//Now the rotation
	ov1 := mol.Coord(orient_atoms[0], 0) //make sure we have the correct versions
	ov2 = mol.Coord(orient_atoms[1], 0)  //same
	orient := ZeroVecs(ov2.NVecs())
	orient.Sub(ov2, ov1)
	//	PDBWrite(mol,"test/2c9v-124centered.pdb")
	Z, _ := NewVecs([]float64{0, 0, 1})
	axis := cross(orient, Z)
	angle := Angle(orient, Z)
	mol.Coords[0] = Rotate(mol.Coords[0], axis, angle)
	if err != nil {
		Te.Error(err)
	}
	PDBWrite("test/2c9v-aligned.pdb", mol.Coords[0], mol, nil)
	fmt.Println("bench1")
}

//TestOldChangeAxis reads the PDB 2c9v.pdb from the test directory, collects
//The CA and CB of residue D124 of the chain A, and rotates the
//whole molecule such as the vector defined by these 2 atoms is
//aligned with the Z axis. The new molecule is written
//as 2c9v_aligned.pdb to the test folder.
func TestOldChangeAxis(Te *testing.T) {
	mol, err := PDBRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	orient_atoms := [2]int{0, 0}
	for index := 0; index < mol.Len(); index++ {
		atom := mol.Atom(index)
		if atom.Chain == "A" && atom.Molid == 124 {
			if atom.Name == "CA" {
				orient_atoms[0] = index
			} else if atom.Name == "CB" {
				orient_atoms[1] = index
			}
		}
	}
	ov1 := mol.Coord(orient_atoms[0], 0)
	ov2 := mol.Coord(orient_atoms[1], 0)
	//now we center the thing in the beta carbon of D124
	mol.Coords[0].SubVec(mol.Coords[0], ov2)
	//Now the rotation
	ov1 = mol.Coord(orient_atoms[0], 0) //make sure we have the correct versions
	ov2 = mol.Coord(orient_atoms[1], 0) //same
	orient := ZeroVecs(ov2.NVecs())
	orient.Sub(ov2, ov1)
	rotation := RotatorToNewZ(orient)
	cr, cc := mol.Coords[0].Dims()
	fmt.Println("rotation: ", rotation, cr, cc) ////////////////////////////////////////////////////////
	mol.Coords[0].Mul(mol.Coords[0], rotation)
	//	fmt.Println(orient_atoms[1], mol.Atom(orient_atoms[1]),mol.Atom(orient_atoms[0]))
	if err != nil {
		Te.Error(err)
	}
	PDBWrite("test/2c9v-old-aligned.pdb", mol.Coords[0], mol, nil)
	fmt.Println("bench2")
}

//Aligns the main plane of a molecule with the XY-plane.
func TestPutInXYPlane(Te *testing.T) {
	mol, err := XYZRead("test/sample_plane.xyz")
	if err != nil {
		Te.Error(err)
	}
	indexes := []int{0, 1, 2, 3, 23, 22, 21, 20, 25, 44, 39, 40, 41, 42, 61, 60, 59, 58, 63, 5}
	some := ZeroVecs(len(indexes))
	some.SomeVecs(mol.Coords[0], indexes)
	//for most rotation things it is good to have the molecule centered on its mean.
	mol.Coords[0], _, _ = MassCentrate(mol.Coords[0], some, nil)
	//The test molecule is not completely planar so we use a subset of atoms that are contained in a plane
	//These are the atoms given in the indexes slice.
	some.SomeVecs(mol.Coords[0], indexes)
	//The strategy is: Take the normal to the plane of the molecule (os molecular subset), and rotate it until it matches the Z-axis
	//This will mean that the plane of the molecule will now match the XY-plane.
	best, err := BestPlane(some, nil)
	if err != nil {
		Te.Error(err)
		panic(err.Error())
	}
	z, _ := NewVecs([]float64{0, 0, 1})
	zero, _ := NewVecs([]float64{0, 0, 0})
	fmt.Println("Best  Plane", best, z)
	axis := ZeroVecs(1)
	axis.Cross(best, z)
	fmt.Println("axis", axis)
	//The main part of the program, where the rotation actually happens. Note that we rotate the whole
	//molecule, not just the planar subset, this is only used to calculate the rotation angle.
	//	fmt.Println("DATA", mol.Coords[0], zero, axis, Angle(best, z))
	mol.Coords[0], err = RotateAbout(mol.Coords[0], zero, axis, Angle(best, z))
	if err != nil {
		Te.Error(err)
	}
	//	fmt.Println("after!", mol.Coords[0], err)
	//Now we write the rotated result.
	XYZWrite("test/Rotated.xyz", mol.Coords[0], mol)
}

func TestDelete(Te *testing.T) {
	mol, err := XYZRead("test/ethanol.xyz")
	if err != nil {
		Te.Error(err)
	}
	mol.Del(4)
	XYZWrite("test/ethanolDel.xyz", mol.Coords[0], mol)

}

func TestWater(Te *testing.T) {
	mol, err := XYZRead("test/sample.xyz")
	if err != nil {
		Te.Error(err)
	}
	for i := 0; i < 6; i++ {
		s := new(Atom)
		if i == 0 || i == 3 {
			s.Symbol = "O"
		} else {
			s.Symbol = "H"
		}
		mol.AppendAtom(s)
	}
	mol.SetCharge(1)
	mol.SetMulti(1)
	c2 := ZeroVecs(mol.Len())
	v := ZeroVecs(6)
	l, _ := mol.Coords[0].Dims()
	fmt.Println(l, mol.Len())
	c2.Stack(mol.Coords[0], v)
	mol.Coords[0] = c2
	c := mol.Coords[0].VecView(43)
	h1 := mol.Coords[0].VecView(42)
	coords := ZeroVecs(mol.Len())
	coords.Copy(mol.Coords[0])
	w1 := MakeWater(c, h1, 2, Deg2Rad*30, true)
	w2 := MakeWater(c, h1, 2, Deg2Rad*-30, false)
	tmp := ZeroVecs(6)
	tmp.Stack(w1, w2)
	fmt.Println("tmp water", w1, w2, tmp, c, h1)
	coords.SetMatrix(mol.Len()-6, 0, tmp)
	XYZWrite("test/WithWater.xyz", coords, mol)
}

func TesstFixPDB(Te *testing.T) {
	mol, err := PDBRead("test/2c9vbroken.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	FixNumbering(mol)
	PDBWrite("test/2c9vfixed.pdb", mol.Coords[0], mol, nil)
}

func TestReduce(Te *testing.T) {
	mol, err := PDBRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	logger, err := os.Create("test/reducereport.log")
	if err != nil {
		Te.Error(err)
	}
	mol2, err := Reduce(mol, mol.Coords[0], 2, logger)
	if err != nil {
		Te.Error(err)
	}
	PDBWrite("test/2c9vHReduce.pdb", mol2.Coords[0], mol2, nil)
}

func TestSuper(Te *testing.T) {
	backbone := []string{"C", "CA", "N"}        //The PDB name of the atoms in the backbone.
	mol1, err := PDBRead("test/2c9v.pdb", true) //true means that we try to read the symbol from the PDB file.
	mol2, err2 := PDBRead("test/1uxm.pdb", true)
	if err != nil || err2 != nil {
		panic("Unable to open input files!")
	}
	mols := []*Molecule{mol1, mol2}
	superlist := make([][]int, 2, 2)
	//We collect the atoms that are part of the backbone.
	for molnumber, mol := range mols {
		for atomindex, atom := range mol.Atoms {
			if isInString(backbone, atom.Name) && atom.Chain == "A" {
				superlist[molnumber] = append(superlist[molnumber], atomindex)
			}
		}
	}
	fmt.Println("superlists!!", len(superlist[0]), len(superlist[1]))
	mol1.Coords[0], err = Super(mol1.Coords[0], mol2.Coords[0], superlist[0], superlist[1])
	if err != nil {
		panic(err.Error())
	}
	newname := "test/2c9v_super.pdb"
	PDBWrite(newname, mol1.Coords[0], mol1, nil)

}
