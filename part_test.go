/// +build part

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
import "os"

//import "time"
import "strings"
import "testing"

//import "os"

//TestMultiXYZ tests that multi-XYZ files are opened and read correctly.
func TestXYZIO(Te *testing.T) {
	mol, err := XYZRead("test/sample.xyz")
	if err != nil {
		fmt.Println("There was an error!")
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
	Z,_ := NewVecs([]float64{0, 0, 1})
	axis := cross(orient, Z)
	angle := AngleInVectors(orient, Z)
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
	rotation := GetSwitchZ(orient)
	cr,cc:=mol.Coords[0].Dims()
	fmt.Println("rotation: ",rotation.Dense, cr, cc) ////////////////////////////////////////////////////////
	mol.Coords[0].Mul(mol.Coords[0],rotation)
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
	}
	z,_ := NewVecs([]float64{0, 0, 1})
	zero,_ := NewVecs([]float64{0, 0, 0})
	fmt.Println("beees", best, z)
	axis := cross(best, z)
	//The main part of the program, where the rotation actually happens. Note that we rotate the whole
	//molecule, not just the planar subset, this is only used to calculate the rotation angle.
	fmt.Println("DATA", mol.Coords[0], zero, axis, AngleInVectors(best, z))
	mol.Coords[0], err = RotateAbout(mol.Coords[0], zero, axis, AngleInVectors(best, z))
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("after!", mol.Coords[0], err)
	//Now we write the rotated result.
	XYZWrite("test/Rotated.xyz", mol.Coords[0], mol)
}

//TestQM tests the QM functionality. It prepares input for ORCA and MOPAC
//In the case of MOPAC it reads a previously prepared output and gets the energy.
func TestQM(Te *testing.T) {
	mol, err := XYZRead("test/sample.xyz")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.Del(mol.Len() - 1)
	mol.SetCharge(1)
	mol.SetMulti(1)
	calc := new(QMCalc)
	calc.SCFTightness = 2 //very demanding
	calc.Optimize = true
	calc.Method = "TPSS"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.HighBasis = "def2-TZVP"
	calc.Grid = 4
	calc.Memory = 1000
	calc.HBAtoms = []int{3, 10, 12}
	calc.HBElements = []string{"Cu", "Zn"}
	calc.CConstraints = []int{0, 10, 20}
	calc.SetDefaults()
	orca := NewOrcaRunner()
	orca.SetnCPU(16) /////////////////////
	atoms, _ := mol.Next(true)
	original_dir, _ := os.Getwd() //will check in a few lines
	if err = os.Chdir("./test"); err != nil {
		Te.Error(err)
	}
	_ = orca.BuildInput(mol, atoms, calc)
	//Now anothertest with HF-3c
	calc.HBAtoms = nil
	calc.HBElements = nil
	calc.RI = false
	calc.Grid = -1
	calc.Dielectric = 0
	calc.Method = "HF-3c"
	orca.SetName("HF3c")
	orca.SetnCPU(8)
	_ = orca.BuildInput(mol, atoms, calc)
	path, _ := os.Getwd()
	//	if err:=orca.Run(false); err!=nil{
	//			Te.Error(err.Error())
	//		}
	fmt.Println(path)
	//Now a MOPAC optimization with the same configuration.
	mopac := NewMopacRunner()
	mopac.BuildInput(mol, atoms, calc)
	mopaccommand := os.Getenv("MOPAC_LICENSE") + "/MOPAC2012.exe"
	mopac.SetCommand(mopaccommand)
	fmt.Println("command", mopaccommand)
	/*	if err := mopac.Run(true); err != nil {
			if strings.Contains(err.Error(), "no such file") {
				fmt.Println("Error", err.Error(), (" Will continue."))
			} else {
				Te.Error(err.Error())
			}
		}
		energy, err := mopac.GetEnergy()
		if err != nil {
			if err.Error() == "Probable problem in calculation" {
				fmt.Println(err.Error())
			} else {
				Te.Error(err)
			}
		}
		geometry, err := mopac.GetGeometry(mol)
		if err != nil {
			if err.Error() == "Probable problem in calculation" {
				fmt.Println(err.Error())
			} else {
				Te.Error(err)
			}
		}
		mol.Coords[0] = geometry
		fmt.Println(energy)
		if err := XYZWrite("mopac.xyz", mol, mol.Coords[0]); err != nil {
			Te.Error(err)
		}
	*/
	//Took away this because it takes too long to run :-)
	/*	if err=orca.Run(true); err!=nil{
		Te.Error(err)
		}
	*/
	if err = os.Chdir(original_dir); err != nil {
		Te.Error(err)
	}
	fmt.Println("end mopac and orca test!")
}

//TestTurbo tests the QM functionality. It prepares input for Turbomole
//Notice that 2 TM inputs cannot be in the same directory. Notice that TMRunner
//supports ECPs
func TesstTurbo(Te *testing.T) {
	mol, err := XYZRead("test/ethanol.xyz")
	os.Chdir("test")
	defer os.Chdir("..")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.Del(mol.Len() - 1)
	mol.SetCharge(0)
	mol.SetMulti(1)
	calc := new(QMCalc)
	calc.SCFConvHelp = 1 //very demanding
	calc.Memory = 1000
	calc.ECP = "ecp-10-mdf"
	calc.ECPElements = []string{"O"}
	calc.Grid = 4
	calc.Optimize = true
	calc.Method = "BP86"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.HighBasis = "def2-TZVP"
	calc.HBElements = []string{"O"}
	calc.RI = true
	calc.Disperssion = "D3"
	calc.CConstraints = []int{0, 3}
	tm := NewTMRunner()
	atoms, _ := mol.Next(true)
	//original_dir, _ := os.Getwd() //will check in a few lines
	//if err = os.Chdir("./test"); err != nil {
	//	Te.Error(err)
	//}
	if err := tm.BuildInput(mol, atoms, calc); err != nil {
		Te.Error(err)
	}
	//os.Chdir(original_dir)
	if err := tm.Run(true); err != nil {
		Te.Error(err)
	}
	energy, err := tm.GetEnergy()
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("energy", energy)
	geo, err := tm.GetGeometry(mol)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("GEO", geo)
	XYZWrite("optiethanol.xyz", geo, mol)
	fmt.Println("end TurboTest!")
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
	coords.Clone(mol.Coords[0])
	w1 := MakeWater(c, h1, 2, Deg2Rad*30, true)
	w2 := MakeWater(c, h1, 2, Deg2Rad*-30, false)
	tmp := ZeroVecs(6)
	tmp.Stack(w1, w2)
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

func TestChemShell(Te *testing.T) {
	mol, err := XYZRead("test/sample.xyz")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.Del(mol.Len() - 1)
	mol.SetCharge(1)
	mol.SetMulti(1)
	calc := new(QMCalc)
	calc.Optimize = true
	calc.Method = "BLYP"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.Grid = 4 //not supported yet, coming sun
	calc.Disperssion = "D3"
	calc.CConstraints = []int{0, 10, 20}
	cs := NewCSRunner()
	atoms, _ := mol.Next(true)
	original_dir, _ := os.Getwd() //will check in a few lines
	if err = os.Chdir("./test"); err != nil {
		Te.Error(err)
	}
	err = cs.BuildInput(mol, atoms, calc)
	qderror_handler(err, Te)
	//now with a PDB
	cs.SetCoordFormat("pdb")
	cs.SetName("gochem_pdb")
	err = cs.BuildInput(mol, atoms, calc)
	qderror_handler(err, Te)
	cs.SetName("gochem_sp")
	calc.Optimize = false
	err = cs.BuildInput(mol, atoms, calc)
	qderror_handler(err, Te)
	if err = os.Chdir(original_dir); err != nil {
		Te.Error(err)
	}
	fmt.Println("end ChemShell test!")
}

func qderror_handler(err error, Te *testing.T) {
	if err != nil {
		if strings.Contains("NonFatal", err.Error()) {
			fmt.Println("Non fatal error: ", err.Error())
		} else {
			Te.Error(err)
		}
	}
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
