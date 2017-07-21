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

import (
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
)

//TestQM tests the QM functionality. It prepares input for ORCA and MOPAC
//In the case of MOPAC it reads a previously prepared output and gets the energy.
func TestQM(Te *testing.T) {
	mol, err := chem.XYZFileRead("../test/sample.xyz")

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
	calc.SetDefaults()
	calc.SCFTightness = 2 //very demanding
	calc.Job = Job{Opti: true}
	//calc.Job.Opti=true
	calc.Method = "TPSS"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.HighBasis = "def2-TZVP"
	calc.Grid = 4
	calc.Memory = 1000
	calc.HBAtoms = []int{3, 10, 12}
	calc.HBElements = []string{"Cu", "Zn"}
	calc.CConstraints = []int{0, 10, 20}
	calc.OldMO = true
	orca := NewOrcaHandle()
	orca.SetnCPU(16) /////////////////////
	atoms := v3.Zeros(mol.Len())
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
		if err := XYZFileWrite("mopac.xyz", mol, mol.Coords[0]); err != nil {
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
func TTestTurbo(Te *testing.T) {
	fmt.Println("Turbomole TEST y wea!")
	mol, err := chem.XYZFileRead("../test/ethanol.xyz")
	original_dir, _ := os.Getwd() //will check in a few lines
	os.Chdir("../test")
	defer os.Chdir(original_dir)
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.SetCharge(0)
	mol.SetMulti(1)
	calc := new(Calc)
	calc.CartesianOpt=true
	calc.SCFConvHelp = 1 //very demanding
	calc.Memory = 1000
	//Not advised
	//	calc.ECP = "ecp-10-mdf"
	//	calc.ECPElements = []string{"O"}
	calc.Grid = 4
	calc.Job = Job{Opti: true}
	calc.Method = "BP86"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.HighBasis = "def2-TZVP"
	calc.HBElements = []string{"C"}
	calc.RI = true
	calc.Dispersion = "D3"
	calc.CConstraints = []int{0, 3}
	tm := NewTMHandle()
	atoms := mol.Coords[0]
	//if err = os.Chdir("./test"); err != nil {
	//	Te.Error(err)
	//}
//	tm.SetDryRun(true) //I don't have TM installed.
	if err := tm.BuildInput(atoms, mol, calc); err != nil {
		Te.Error(err)
	}
	/*

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
		chem.XYZFileWrite("optiethanol.xyz", geo, mol)
		fmt.Println("end TurboTest!")
	*/
	//	os.Chdir(original_dir)
}

/*
func TestFermions(Te *testing.T) {
	mol, err := chem.XYZFileRead("../test/ethanol.xyz")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	mol.SetCharge(0)
	mol.SetMulti(1)
	calc := new(Calc)
	calc.Job=Job{Opti:true}
	calc.Method = "BLYP"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.Grid = 4
	calc.Dispersion = "D3"
	calc.CConstraints = []int{0, 10, 20}
	cs := NewFermionsHandle()
	cs.SetName("gochemF")
	atoms := mol.Coords[0]
	original_dir, _ := os.Getwd() //will check in a few lines
	if err = os.Chdir("../test"); err != nil {
		Te.Error(err)
	}
	err = cs.BuildInput(atoms, mol, calc)
	defer os.Chdir(original_dir)
	E, err := cs.Energy()
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("Final energy:", E, "kcal/mol")
	//	ngeo,err:=cs.OptimizedGeometry(mol)
	//	if err!=nil{
	//		fmt.Println("Error with the geometry?: ", err.Error())
	//	}
	//	chem.XYZFileWrite("LastGeoFermions.xyz",ngeo,mol)
	fmt.Println("Passed FermiONs++ test!")
}
*/

func qderror_handler(err error, Te *testing.T) {
	if err != nil {
		if strings.Contains("NonFatal", err.Error()) {
			fmt.Println("Non fatal error: ", err.Error())
		} else {
			Te.Error(err)
		}
	}
}

func TestNWChem(Te *testing.T) {
	mol, err := chem.XYZFileRead("../test/ethanol.xyz")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	fmt.Println(mol.Coords[0], len(mol.Coords), "Quiere quedar leyenda, compadre?", err)
	mol.SetCharge(0)
	mol.SetMulti(1)
	calc := new(Calc)
	calc.SCFTightness = 1 //quite tight
	calc.SCFConvHelp = 1
	calc.Job = Job{Opti: true}
	calc.Method = "TPSS"
	calc.Dielectric = 4
	calc.Basis = "def2-SVP"
	calc.HighBasis = "def2-TZVP"
	calc.Grid = 4
	calc.Memory = 1000
	calc.HBAtoms = []int{2}
	calc.HBElements = []string{"O"}
	calc.CConstraints = []int{1}
	calc.SetDefaults()
	nw := NewNWChemHandle()
	orca := NewOrcaHandle()
	nw.SetName("gochem")
	orca.SetName("gochemII")
	atoms := v3.Zeros(mol.Len())
	mol.Next(atoms)
	if err = os.Chdir("../test"); err != nil {
		Te.Error(err)
	}
	err = nw.BuildInput(atoms, mol, calc)
	if err != nil {
		Te.Error(err)
	}
	_ = orca.BuildInput(atoms, mol, calc)
	//The files are already in ./test.
	os.Chdir("../test")
	defer os.Chdir("../qm")
	energy, err := nw.Energy()
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("NWChem Energy: ", energy)
	newg, err := nw.OptimizedGeometry(mol)
	if err != nil {
		Te.Error(err)
	}
	chem.XYZFileWrite("optiNW.xyz", newg, mol)

}


func TestXtb (Te *testing.T) {
	mol, err := chem.XYZFileRead("../test/ethanol.xyz")
	if err != nil {
		Te.Error(err)
	}
	if err := mol.Corrupted(); err != nil {
		Te.Error(err)
	}
	fmt.Println(mol.Coords[0], len(mol.Coords), "Quiere quedar XTB leyenda, compadre?", err)
	mol.SetCharge(0)
	mol.SetMulti(1)
	calc := new(Calc)
	calc.Job = Job{Opti: true}
	//no support for constraints yet
	calc.Method = "" //we only use xtb here soooo
	calc.Dielectric = 4
	xtb := NewXTBHandle()
	xtb.SetName("XTBgochem")
	atoms := v3.Zeros(mol.Len())
	mol.Next(atoms)
	if err = os.Chdir("../test"); err != nil {
		Te.Error(err)
	}
	err = xtb.BuildInput(atoms, mol, calc)
	if err != nil {
		Te.Error(err)
	}
	if err := xtb.Run(true); err != nil {
		Te.Error(err)
	}
	os.Chdir("../test")
	defer os.Chdir("../qm")
	energy, err := xtb.Energy()
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("XTB Energy: ", energy)
	newg, err := xtb.OptimizedGeometry(mol)
	if err != nil {
		Te.Error(err)
	}
	chem.XYZFileWrite("optiXTB.xyz", newg, mol)

}
