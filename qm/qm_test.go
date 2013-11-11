/*
 * qm_test.go, part of gochem.
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

package qm

import(
	"testing"
	"github.com/rmera/gochem"
	"fmt"
	"os"
	"strings"
)


//TestQM tests the QM functionality. It prepares input for ORCA and MOPAC
//In the case of MOPAC it reads a previously prepared output and gets the energy.
func TestQM(Te *testing.T) {
	mol, err := chem.XYZRead("../test/sample.xyz")

	fmt.Println(mol.Coords[0], len(mol.Coords), "LOS JUIMOS CTM", err)
	if err != nil {
		Te.Error(err)

	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.Del(mol.Len() - 1)
	mol.SetCharge(1)
	mol.SetMulti(1)
	calc := new(Calc)
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
	orca := NewOrcaHandle()
	orca.SetnCPU(16) /////////////////////
	atoms := chem.ZeroVecs(mol.Len())
	mol.Next(atoms)
	original_dir, _ := os.Getwd() //will check in a few lines
	if err = os.Chdir("../test"); err != nil {
		Te.Error(err)
	}
	_ = orca.BuildInput(atoms, mol, calc)
	//Now anothertest with HF-3c
	calc.HBAtoms = nil
	calc.HBElements = nil
	calc.RI = false
	calc.Grid = -1
	calc.Dielectric = 0
	calc.Method = "HF-3c"
	orca.SetName("HF3c")
	orca.SetnCPU(8)
	//	fmt.Println(mol.Coords[0], "vieja")
	_ = orca.BuildInput(atoms, mol, calc)
	path, _ := os.Getwd()
	//	if err:=orca.Run(false); err!=nil{
	//			Te.Error(err.Error())
	//		}
	fmt.Println(path)
	//Now a MOPAC optimization with the same configuration.
	mopac := NewMopacHandle()
	mopac.BuildInput(atoms, mol, calc)
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
		energy, err := mopac.Energy()
		if err != nil {
			if err.Error() == "Probable problem in calculation" {
				fmt.Println(err.Error())
			} else {
				Te.Error(err)
			}
		}
		geometry, err := mopac.OptimizedGeometry(mol)
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
//Notice that 2 TM inputs cannot be in the same directory. Notice that TMHandle
//supports ECPs
func TesstTurbo(Te *testing.T) {
	mol, err := chem.XYZRead("../test/ethanol.xyz")
	os.Chdir("test")
	defer os.Chdir("..")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.SetCharge(0)
	mol.SetMulti(1)
	calc := new(Calc)
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
	tm := NewTMHandle()
	atoms := mol.Coords[0]
	//original_dir, _ := os.Getwd() //will check in a few lines
	//if err = os.Chdir("./test"); err != nil {
	//	Te.Error(err)
	//}
	if err := tm.BuildInput(atoms, mol, calc); err != nil {
		Te.Error(err)
	}
	//os.Chdir(original_dir)
	if err := tm.Run(true); err != nil {
		Te.Error(err)
	}
	energy, err := tm.Energy()
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("energy", energy)
	geo, err := tm.OptimizedGeometry(mol)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("GEO", geo)
	chem.XYZWrite("optiethanol.xyz", geo, mol)
	fmt.Println("end TurboTest!")
}



func TestChemShell(Te *testing.T) {
	mol, err := chem.XYZRead("../test/sample.xyz")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.Del(mol.Len() - 1)
	mol.SetCharge(1)
	mol.SetMulti(1)
	calc := new(Calc)
	calc.Optimize = true
	calc.Method = "BLYP"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.Grid = 4 //not supported yet, coming sun
	calc.Disperssion = "D3"
	calc.CConstraints = []int{0, 10, 20}
	cs := NewCSHandle()
	atoms := mol.Coords[0]
	original_dir, _ := os.Getwd() //will check in a few lines
	if err = os.Chdir("../test"); err != nil {
		Te.Error(err)
	}
	err = cs.BuildInput(atoms, mol, calc)
	qderror_handler(err, Te)
	//now with a PDB
	cs.SetCoordFormat("pdb")
	cs.SetName("gochem_pdb")
	err = cs.BuildInput(atoms, mol, calc)
	qderror_handler(err, Te)
	cs.SetName("gochem_sp")
	calc.Optimize = false
	err = cs.BuildInput(atoms, mol, calc)
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

