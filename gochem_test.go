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

import (
	"fmt"
	"os"
	"runtime"
	"testing"

	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

//import "runtime"

func TestGROIO(Te *testing.T) {
	mol, err := GroFileRead("test/test.gro")
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("NO WEI", mol.Len(), len(mol.Coords))
	fmt.Println(mol.Atom(3), mol.Coords[0].VecView(3))
	err = PDBFileWrite("test/testgro.pdb", mol.Coords[0], mol, nil)
	if err != nil {
		Te.Error(err)
	}
	err = GroFileWrite("test/testgro.gro", mol.Coords, mol)
	if err != nil {
		Te.Error(err)
	}

}

//TestMultiXYZ tests that multi-XYZ files are opened and read correctly.
func TestXYZIO(Te *testing.T) {
	mol, err := XYZFileRead("test/sample.xyz")
	if err != nil {
		fmt.Println("There was an error!", err.Error())
		Te.Error(err)
	}
	fmt.Println("XYZ read!")
	XYZFileWrite("test/sampleFirst.xyz", mol.Coords[0], mol)
}

func TestPDBIO(Te *testing.T) {
	mol, err := PDBFileRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("NO WEI")
	err = PDBFileWrite("test/2c9vIO.pdb", mol.Coords[0], mol, mol.Bfactors[0])
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
	//runtime.GOMAXPROCS(2) ///////////////////////////
	mol, err := PDBFileRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	PDBFileWrite("test/2c9v-Readtest.pdb", mol.Coords[0], mol, nil)
	//The selection thing
	orient_atoms := [2]int{0, 0}
	for index := 0; index < mol.Len(); index++ {
		atom := mol.Atom(index)
		if atom.Chain == "A" && atom.MolID == 124 {
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
	PDBFileWrite("test/2c9v-translated.pdb", mol.Coords[0], mol, nil)
	//Now the rotation
	ov1 := mol.Coord(orient_atoms[0], 0) //make sure we have the correct versions
	ov2 = mol.Coord(orient_atoms[1], 0)  //same
	orient := v3.Zeros(ov2.NVecs())
	orient.Sub(ov2, ov1)
	//	PDBWrite(mol,"test/2c9v-124centered.pdb")
	Z, _ := v3.NewMatrix([]float64{0, 0, 1})
	axis := cross(orient, Z)
	angle := Angle(orient, Z)
	oldcoords := v3.Zeros(mol.Coords[0].NVecs())
	oldcoords.Copy(mol.Coords[0])
	mol.Coords[0] = Rotate(oldcoords, mol.Coords[0], axis, angle)
	if err != nil {
		Te.Error(err)
	}
	PDBFileWrite("test/2c9v-aligned.pdb", mol.Coords[0], mol, nil)
	fmt.Println("bench1")
}

func TestMolidNameChain2Index(Te *testing.T) {
	//runtime.GOMAXPROCS(2) ///////////////////////////
	mol, err := PDBFileRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	index, err := MolIDNameChain2Index(mol, 46, "ND1", "A")
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("this should print the index and the Atom object for H46, chain A: ", index, mol.Atom(index))

}

//TestOldChangeAxis reads the PDB 2c9v.pdb from the test directory, collects
//The CA and CB of residue D124 of the chain A, and rotates the
//whole molecule such as the vector defined by these 2 atoms is
//aligned with the Z axis. The new molecule is written
//as 2c9v_aligned.pdb to the test folder.
func TestOldChangeAxis(Te *testing.T) {
	viej, _ := os.Open("test/2c9v.pdb")
	mol, err := PDBRead(viej, true)
	if err != nil {
		Te.Error(err)
	}
	orient_atoms := [2]int{0, 0}
	for index := 0; index < mol.Len(); index++ {
		atom := mol.Atom(index)
		if atom.Chain == "A" && atom.MolID == 124 {
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
	orient := v3.Zeros(ov2.NVecs())
	orient.Sub(ov2, ov1)
	rotation := RotatorToNewZ(orient)
	cr, cc := mol.Coords[0].Dims()
	fmt.Println("rotation: ", rotation, cr, cc) ////////////////////////////////////////////////////////
	mol.Coords[0].Mul(mol.Coords[0], rotation)
	//	fmt.Println(orient_atoms[1], mol.Atom(orient_atoms[1]),mol.Atom(orient_atoms[0]))
	if err != nil {
		Te.Error(err)
	}
	PDBFileWrite("test/2c9v-old-aligned.pdb", mol.Coords[0], mol, nil)
	fmt.Println("bench2")
}

//Aligns the main plane of a molecule with the XY-plane.
//Here XYZRead and XYZWrite are tested
func TestPutInXYPlane(Te *testing.T) {
	myxyz, _ := os.Open("test/sample_plane.xyz")
	mol, err := XYZRead(myxyz)
	if err != nil {
		Te.Error(err)
	}
	indexes := []int{0, 1, 2, 3, 23, 22, 21, 20, 25, 44, 39, 40, 41, 42, 61, 60, 59, 58, 63, 5}
	some := v3.Zeros(len(indexes))
	some.SomeVecs(mol.Coords[0], indexes)
	//for most rotation things it is good to have the molecule centered on its mean.
	mol.Coords[0], _, _ = MassCenter(mol.Coords[0], some, nil)
	//The test molecule is not completely planar so we use a subset of atoms that are contained in a plane
	//These are the atoms given in the indexes slice.
	some.SomeVecs(mol.Coords[0], indexes)
	//The strategy is: Take the normal to the plane of the molecule (os molecular subset), and rotate it until it matches the Z-axis
	//This will mean that the plane of the molecule will now match the XY-plane.
	best, err := BestPlane(some, nil)
	if err != nil {
		err2 := err.(Error)
		fmt.Println(err2.Decorate(""))
		Te.Errorf(err.Error())
		//		panic(err.Error())
	}
	z, _ := v3.NewMatrix([]float64{0, 0, 1})
	zero, _ := v3.NewMatrix([]float64{0, 0, 0})
	fmt.Println("Best  Plane", best, z)
	axis := v3.Zeros(1)
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
	outxyz, _ := os.Create("test/Rotated.xyz") //This is the XYZWrite written file
	XYZWrite(outxyz, mol.Coords[0], mol)
}

func TestDelete(Te *testing.T) {
	mol, err := XYZFileRead("test/ethanol.xyz")
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("Calling with 8")
	mol.Del(8)
	XYZFileWrite("test/ethanolDel8.xyz", mol.Coords[0], mol)
	mol2, err := XYZFileRead("test/ethanol.xyz")
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("Calling with 4")
	mol2.Del(4)
	XYZFileWrite("test/ethanolDel4.xyz", mol2.Coords[0], mol2)

}

func TestWater(Te *testing.T) {
	//	runtime.GOMAXPROCS(2) ///////////////////////////
	mol, err := XYZFileRead("test/sample.xyz")
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
	c2 := v3.Zeros(mol.Len())
	v := v3.Zeros(6)
	l, _ := mol.Coords[0].Dims()
	fmt.Println(l, mol.Len())
	c2.Stack(mol.Coords[0], v)
	mol.Coords[0] = c2
	c := mol.Coords[0].VecView(43)
	h1 := mol.Coords[0].VecView(42)
	coords := v3.Zeros(mol.Len())
	coords.Copy(mol.Coords[0])
	w1 := MakeWater(c, h1, 2, Deg2Rad*30, true)
	w2 := MakeWater(c, h1, 2, Deg2Rad*-30, false)
	tmp := v3.Zeros(6)
	tmp.Stack(w1, w2)
	fmt.Println("tmp water", w1, w2, tmp, c, h1)
	coords.SetMatrix(mol.Len()-6, 0, tmp)
	XYZFileWrite("test/WithWater.xyz", coords, mol)
	fmt.Println("Done TestWater")
}

func TesstFixPDB(Te *testing.T) {
	mol, err := PDBFileRead("test/2c9vbroken.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	FixNumbering(mol)
	PDBFileWrite("test/2c9vfixed.pdb", mol.Coords[0], mol, nil)
	fmt.Println("DoneTestFixPDB")
}

//will fail if reduce is not installed!
func TTTestReduce(Te *testing.T) { //silenced
	fmt.Println("Start TestReduce")
	mol, err := PDBFileRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	logger, err := os.Create("test/reducereport.log")
	if err != nil {
		Te.Error(err)
	}
	defer logger.Close()
	mol2, err := Reduce(mol, mol.Coords[0], 2, logger, "")
	if err != nil {
		Te.Error(err)
	}
	PDBFileWrite("test/2c9vHReduce.pdb", mol2.Coords[0], mol2, nil)
	fmt.Println("END TestReduce")
}

func TTestShape(Te *testing.T) {
	myhandle, _ := os.Open("test/2c9v.pdb")
	mol1, err := PDBRead(myhandle, true) //true means that we try to read the symbol from the PDB file.
	masses, err := mol1.Masses()
	if err != nil {
		Te.Error(err)
	}
	moment, err := MomentTensor(mol1.Coords[0], masses)
	if err != nil {
		Te.Error(err)
	}
	rhos, err := Rhos(moment)
	if err != nil {
		Te.Error(err)
	}
	linear, circular, err := RhoShapeIndexes(rhos)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("liner,circular distortion:", linear, circular)
	lin2, circ2, err := EasyShape(mol1.Coords[0], -1, mol1)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("Easier way linear,circular:", lin2, circ2)
	mol2, _ := XYZFileRead("test/sample_plane.xyz")
	lin3, circ3, err := EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("sample_plane.xyz shape indicators; linear,circular:", lin3, circ3)
	//now the shapetests batterty!
	mol2, _ = XYZFileRead("test/shapetests/porphyrin.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("porphyrin.xyz shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/2-mesoporphyrin.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("2-mesoporphyrin.xyz shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/4-mesoporphyrin.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("4-mesoporphyrin.xyz shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/heptane.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("heptane.xyz shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/decane.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("decane.xyz shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/phenantrene.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("phenantrene.xyz shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/methylphenantrene.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("methylphenantrene shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/tbutylphenantrene.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("tbutylphenantrene shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/fullerene20.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], -1, mol2)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("fullerene20.xyz shape indicators; linear,circular:", lin3, circ3)
	mol2, _ = XYZFileRead("test/shapetests/fullerene60.xyz")
	lin3, circ3, err = EasyShape(mol2.Coords[0], 0.0001, mol2) //maybe it's too symmetrical for the default epsilon?
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("fullerene60.xyz shape indicators; linear,circular:", lin3, circ3)

}

//Here PDBRead and PDBWrite are tested
func TestSuper(Te *testing.T) {
	backbone := []string{"CA", "C", "N"} //The PDB name of the atoms in the backbone.
	myhandle, _ := os.Open("test/2c9v.pdb")
	mol1, err := PDBRead(myhandle, true) //true means that we try to read the symbol from the PDB file.
	mol2, err2 := PDBFileRead("test/1uxm.pdb", true)
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
				//	fmt.Println(atom)
			}
		}
	}
	fmt.Println("superlists!!", len(superlist[0]), len(superlist[1]))
	mol1.Coords[0], err = Super(mol1.Coords[0], mol2.Coords[0], superlist[0], superlist[1])
	rmsd1, _ := rMSD(mol1.Coords[0], mol2.Coords[0], superlist[0], superlist[1])
	rmsd2, _ := RMSD(mol1.Coords[0], mol2.Coords[0], superlist[0], superlist[1])
	fmt.Println("RMSDs for proteins!", rmsd2, rmsd1)
	fmt.Println("Atoms superimposed:", len(superlist[0]))
	if err != nil {
		panic(err.Error())
	}
	newname := "test/2c9v_super.pdb" //This is the PDBWrite written file
	pdbout, _ := os.Create(newname)
	PDBWrite(pdbout, mol1.Coords[0], mol1, nil)
	//Now for a full molecule
	ptest, _ := XYZFileRead("test/Rotated.xyz")
	ptempla, _ := XYZFileRead("test/sample_plane.xyz")
	newp, err := Super(ptest.Coords[0], ptempla.Coords[0], nil, nil)
	rmsd2, _ = RMSD(newp, ptempla.Coords[0])
	rmsd3, _ := RMSD(newp, ptempla.Coords[0], nil, nil)
	rmsd1, _ = rMSD(newp, ptempla.Coords[0], nil, nil)
	fmt.Println("RMSD mol (should be 0):", rmsd1, rmsd2, rmsd3)
	if err != nil {
		panic(err.Error())
	}
	XYZFileWrite("test/SuperPlane.xyz", newp, ptest)

}

func TestRotateBz(Te *testing.T) {
	runtime.GOMAXPROCS(2)
	fmt.Println("Here we go!")
	mol, err := XYZFileRead("test/BZ.xyz")
	if err != nil {
		panic(err.Error())
	}
	carbonIn := []int{}
	bzIn := []int{}
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if at.Symbol == "C" {
			bzIn = append(bzIn, i)
			carbonIn = append(carbonIn, i)
		} else if at.Symbol == "H" {
			bzIn = append(bzIn, i)
		}
	}
	coordsI := mol.Coords[0]
	carbons := v3.Zeros(6)
	bz := v3.Zeros(12)
	carbons.SomeVecs(coordsI, carbonIn)
	coords := v3.Zeros(mol.Len())
	coords, _, _ = MassCenter(coordsI, carbons, nil)
	bz.SomeVecs(coords, bzIn)
	carbons.SomeVecs(coords, carbonIn)
	planevec, err := BestPlane(carbons, nil)
	if err != nil {
		if e, ok := err.(Error); ok {
			fmt.Println("DEcoration:", e.Decorate(""))
		}
		Te.Errorf(err.Error())
	}
	basename := "BZ"
	newcoords := v3.Zeros(mol.Len())
	origin := v3.Zeros(1)
	bzcopy := v3.Zeros(12)
	bzcopy2 := v3.Zeros(12) //testing
	rot := v3.Zeros(12)
	rot3 := v3.Zeros(12)
	for _, angle := range []float64{0, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 180} {
		bzcopy.Copy(bz)
		bzcopy2.Copy(bz) //testing
		rot = Rotate(bzcopy, rot, planevec, Deg2Rad*angle)
		rot3 = RotateSer(bzcopy, rot, planevec, Deg2Rad*angle)
		rot2, _ := EulerRotateAbout(bzcopy2, origin, planevec, Deg2Rad*angle) //should be the same as the previous
		if !mat.EqualApprox(rot, rot2, 0.01) {
			Te.Errorf("Rotors Rotate and EulerRotate not equal for angle %3.2f", angle)
		} else if !mat.EqualApprox(rot2, rot3, 0.01) {
			Te.Errorf("Rotors RotateSer and EulerRotate not equal for angle %3.2f", angle)

		} else {
			fmt.Println("Rotors EQUAL for angle", angle)

		}
		fmt.Println("rot", rot, "rot2", rot2)
		newcoords.Copy(coords)
		newcoords.SetVecs(rot, bzIn)
		//test
		//	tempcoords.Stack(planevec,origin)
		//	testxyz.Stack(newcoords,tempcoords)
		//end
		XYZFileWrite(fmt.Sprintf("test/%s-%3.1f.xyz", basename, angle), newcoords, mol)

	}
	//	fmt.Println(mol, planevec)
}

func TestProjectionAndAntiProjection(Te *testing.T) {
	A := v3.Zeros(1)
	A.Set(0, 0, 2.0)
	B, _ := v3.NewMatrix([]float64{1, 1, 0})
	C := AntiProjection(A, B)
	D := Projection(B, A)
	fmt.Println("Projection of B on A (D)", D)
	fmt.Println("Anti-projection of A on B (C):", C)
	fmt.Println("Norm of C: ", C.Norm(0), " Norm of A,B: ", A.Norm(0), B.Norm(0), "Norm of D:", D.Norm(0))
}
