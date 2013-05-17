// +build part

/*
 * part_test.go
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


//import "github.com/skelterjohn/go.matrix"
import "fmt"
import "testing"
//import "os"

	
//TestMultiXyz tests that multi-XYZ files are opened and read correctly.
func TestXyzIO(Te *testing.T) {
	mol, err := XyzRead("test/sample.xyz")
	if err != nil {
		fmt.Println("There was an error!")
		Te.Error(err)
	}
	fmt.Println("XYZ read!")
	XyzWrite("test/sampleFirst.xyz", mol, mol.Coords[0])
}

func TestPDBIO(Te *testing.T) {
	mol, err := PdbRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	err=PdbWrite("test/2c9vIO.pdb", mol, mol.Coords[0],mol.Bfactors[0])
	if err != nil {
		Te.Error(err)
	}
	//for the 3 residue  I should get -131.99, 152.49.
}



//TestOldChangeAxis reads the PDB 2c9v.pdb from the test directory, collects 
//The CA and CB of residue D124 of the chain A, and rotates the 
//whole molecule such as the vector defined by these 2 atoms is 
//aligned with the Z axis. The new molecule is written
//as 2c9v_aligned.pdb to the test folder.
func TestOldChangeAxis(Te *testing.T) {
	mol, err := PdbRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	orient_atoms := [2]int{0, 0}
	for index := 0; index < mol.Len(); index++ {
		atom := mol.Atom(index)
		if atom.Chain == 'A' && atom.Molid == 124 {
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
	mol.Coords[0].SubRow(mol.Coords[0], ov2)
	//Now the rotation
	ov1 = mol.Coord(orient_atoms[0], 0) //make sure we have the correct versions
	ov2 = mol.Coord(orient_atoms[1], 0) //same
	orient := gnClone(ov2)
	orient.Sub(orient,ov1)
	rotation := GetSwitchZ(orient)
	//	fmt.Println("rotation: ",rotation)
	mol.Coords[0] = gnMul(mol.Coords[0], rotation)
	//	fmt.Println(orient_atoms[1], mol.Atom(orient_atoms[1]),mol.Atom(orient_atoms[0]))
	if err != nil {
		Te.Error(err)
	}
	PdbWrite("test/2c9v-old-aligned.pdb", mol, mol.Coords[0])
	fmt.Println("bench2")
}


//Aligns the main plane of a molecule with the XY-plane.
func TestPutInXYPlane(Te *testing.T) {
	mol, err := XyzRead("test/sample_plane.xyz")
	if err != nil {
		Te.Error(err)
	}
	indexes := []int{0, 1, 2, 3, 23, 22, 21, 20, 25, 44, 39, 40, 41, 42, 61, 60, 59, 58, 63, 5}
	some:=gnZeros(len(indexes),3)
	some.SomeRows(mol.Coords[0], indexes)
	//for most rotation things it is good to have the molecule centered on its mean.
	mol.Coords[0], _, _ = MassCentrate(mol.Coords[0], some, nil)
	//The test molecule is not completely planar so we use a subset of atoms that are contained in a plane
	//These are the atoms given in the indexes slice.
	some.SomeRows(mol.Coords[0], indexes)
	//The strategy is: Take the normal to the plane of the molecule (os molecular subset), and rotate it until it matches the Z-axis
	//This will mean that the plane of the molecule will now match the XY-plane.
	best, err := BestPlane(nil, some)
	if err != nil {
		Te.Error(err)
	}
	z := NewCoord([]float64{0, 0, 1}, 1, 3)
	zero := NewCoord([]float64{0, 0, 0}, 1, 3)
	axis, err := Cross3D(best, z)
	if err != nil {
		Te.Error(err)
	}
	//The main part of the program, where the rotation actually happens. Note that we rotate the whole
	//molecule, not just the planar subset, this is only used to calculate the rotation angle.
	mol.Coords[0], err = RotateAbout(mol.Coords[0], zero, axis, AngleInVectors(best, z))
	if err != nil {
		Te.Error(err)
	}
	//Now we write the rotated result.
	XyzWrite("test/Rotated.xyz", mol, mol.Coords[0])
}
