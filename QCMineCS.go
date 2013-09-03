
/*
 * qcminecs.go, part of gochem.
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

//Interface to run QM jobs in QCMine/Chemshell

package chem

import "os"
import "io"
import "strings"

//import "strconv"
import "bufio"
import "fmt"
import "os/exec"


type QCCSRunner struct {
	defmethod   string
	defbasis    string
	coordformat string
	previousMO  string
	command     string
	inputname   string
	gimic       bool
}

//Creates and initialized a new instance of QCCSRuner, with values set
//to its defaults.
func MakeQCCSRunner() *QCCSRunner {
	run := new(QCCSRunner)
	run.SetDefaults()
	return run
}

//QCCSRunner methods

//Just to satisfy the interface. It does nothing
func (O *QCCSRunner) SetnCPU(cpu int) {
	//It does nothing! :-D
}

//This set the name of the subdirectory, in the current directory
//where the calculation will be ran
func (O *QCCSRunner) SetName(name string) {
	O.inputname = name

}

func (O *QCCSRunner) SetCoordFormat(format string){
	f,ok :=chemShellFormats[format]
	if !ok{
		f="xyz"
	}
	O.coordformat=f
}

var chemShellFormats  = map[string]string{
	"XYZ": "xyz",
	"xyz": "xyz",
	"PDB": "pdb",
	"pdb": "pdb",
	
}


//SetCommand sets the command to run the ChemShell/QCMine calculation.
func (O *QCCSRunner) SetCommand(name string) {
	//Does nothing again
}

//Sets some defaults for QCCSRunner. default is an optimization at
//  PBE0-D3-gCP / def2-SVP
func (O *QCCSRunner) SetDefaults() {
	O.defmethod = "pbe0-d"
	O.defbasis = "def2-SVP"
	O.inputname= "qcmine"
	O.coordformat="xyz"
}


func (O *QCCSRunner) BuildInput(atoms Ref, coords *CoordMatrix, Q *QMCalc) error {
	basis,ok:=qcMineBasis[Q.Basis]
	if !ok {
		basis=O.defbasis
	}
	method,ok:=qcMineMethods[Q.Method]
	if !ok{
		method=O.defmethod
	}
	disp,ok:=qcMineDisp[Q.Disperssion]
	if !ok {
		disp="d"
	}
    if !strings.HasSuffix(method,"-d") && disp!=""{
		method=strings.Join([]string{method,disp},"-")
	}
	if Q.BSSE=="gCP"{
		method=strings.Join([]string{method,"gcp"},"-")
	}
	coordline:=""
	if O.coordformat=="pdb"{
		PDBWrite("coord.pdb",atoms,coords)
		coordline="read_pdb coord.pdb coords=coord.crd\n"
	}else{
		XYZWrite("coord.xyz",atoms,coords)
		coordline="read_xyz coord.xyz coords=coord.crd\n"
		
	}

}




var qcMineDisp = map[string]string{
	"D3": "d",
}


var qcMineMethods = map[string]string{
	"HF":     "hf",
	"hf":     "hf",
	"b3lyp":  "b3-lyp",
	"B3LYP":  "b3-lyp",
	"b3-lyp": "b3-lyp",
	"PBE0":    "pbe0",
	"pbe0":    "pbe0",
	"BP86":   "bp86",
	"bp86":    "bp86",
	"revpbe":  "revPBE-d",
	"revPBE":   "revPBE-d",
}

var qcMineBasis = map[string]string{
	"def2-TZVP": "tzp",
	"tzp": "tzp",
	"def2-SVP": "svp",
	"svp": "svp",
}



