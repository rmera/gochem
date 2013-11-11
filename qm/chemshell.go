/*
 * chemshell.go, part of gochem.
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
 * Gochem was started at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland. It is currently developed at AK Ochsenfeld, LMU-Munich,
 * Germany.
 *
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

//Interface to run QM jobs in QCMine/Chemshell

package qm

import (
	"fmt"
	"os"
	"strings"
	"github.com/rmera/gochem"
)

type CSHandle struct {
	program     string
	defmethod   string
	defbasis    string
	coordformat string
	previousMO  string
	command     string
	inputname   string
	link        bool
	gimic       bool
	mpi			bool
}

//Creates and initialized a new instance of QCCSRuner, with values set
//to its defaults.
func NewCSHandle() *CSHandle {
	run := new(CSHandle)
	run.SetDefaults()
	return run
}

//QCCSHandle methods

//Just to satisfy the interface. It does nothing
func (O *CSHandle) SetnCPU(cpu int) {
	//It does nothing! :-D
}

func (O *CSHandle) SetMPI(mpi bool){
	O.mpi=mpi
}

//This set the name of the subdirectory, in the current directory
//where the calculation will be ran
func (O *CSHandle) SetName(name string) {
	O.inputname = name

}

func (O *CSHandle) SetCoordFormat(format string) {
	f, ok := chemShellFormats[format]
	if !ok {
		f = "xyz"
	}
	O.coordformat = f
}

var chemShellFormats = map[string]string{
	"XYZ": "xyz",
	"xyz": "xyz",
	"PDB": "pdb",
	"pdb": "pdb",
}

//SetCommand sets the command to run the ChemShell/QCMine calculation.
func (O *CSHandle) SetCommand(name string) {
	//Does nothing again
}

//Sets some defaults for QCCSHandle. default is an optimization at
//  PBE0-D3-gCP / def2-SVP
func (O *CSHandle) SetDefaults() {
	O.defmethod = "pbe0-d"
	O.defbasis = "def2-SVP"
	O.link = true
	O.inputname = "gochem"
	O.coordformat = "xyz"
	O.program = "qcmine"
}

//BuildInput builds a ChemShell input (at this point only for pure QM calculations with QCMine). Returns error on failure.
func (O *CSHandle) BuildInput(coords *chem.VecMatrix, atoms chem.ReadRef, Q *Calc) error {
	var nonfatal error
	if atoms.Multi() != 1 {
		return fmt.Errorf("Only closed shell supported for ChemShell")
	}
	if O.program != "qcmine" {
		nonfatal = fmt.Errorf("NonFatal: Unavailable program requested. Only QCMine is supported at this point. Will switch to QCMine")
		O.program = "qcmine"
	}
	basis, ok := chemShellBasis[Q.Basis]
	if !ok {
		nonfatal = fmt.Errorf("NonFatal: Unavailable basis set requested, using default")
		basis = O.defbasis
	}
	method, ok := chemShellMethods[Q.Method]
	if !ok {
		nonfatal = fmt.Errorf("NonFatal: Unavailable method requested: %s, using default: %s",Q.Method,O.defmethod)
		method = O.defmethod
	}
	disp, ok := qcMineDisp[Q.Disperssion]
	if !ok {
		nonfatal = fmt.Errorf("NonFatal: Unavailable disperssion correction requested, using default")
		disp = "d"
	}
	if !strings.HasSuffix(method, "-d") && disp != "" {
		method = strings.Join([]string{method, disp}, "-")
	}
	if Q.BSSE == "gCP" {
		method = strings.Join([]string{method, "gcp"}, "-")
	}
	coordline := ""
	if O.coordformat == "pdb" {
		chem.PDBWrite(fmt.Sprintf("%s.pdb", O.inputname), coords, atoms, nil)
		coordline = fmt.Sprintf("set residues [ pdb_to_res %s.pdb ] \nread_pdb file=%s.pdb coords=%s.crd\n", O.inputname, O.inputname, O.inputname)

	} else {
		chem.XYZWrite(fmt.Sprintf("%s.xyz", O.inputname), coords, atoms)
		coordline = fmt.Sprintf("read_xyz %s.xyz coords=%s.crd \nset residues [ res_selectall coords=%s.crd ]\n", O.inputname, O.inputname, O.inputname)
	}
	frozen := ""
	optline := ""
	writeline := ""
	if Q.Optimize {
		optline = "coordinates=hdlc \\\n     "
		if len(Q.CConstraints) > 0 {
			frozen = "set frozen [ list "
			for _, v := range Q.CConstraints {
				frozen = fmt.Sprintf("%s %d", frozen, v)
			}
			frozen = fmt.Sprintf("%s]\n", frozen)
			optline = fmt.Sprintf("%s frozen= $frozen \\\n     ", optline)
		}
		optline = fmt.Sprintf("%s residues= $residues \\\n     maxcycle=300 \\\n     result=%s.crd \\\n     ", optline, O.inputname)
		writeline = fmt.Sprintf("write_xyz file=%s.xyz coords=%s.crd\n", O.inputname, O.inputname)
	}
	//we start actually writing the input
	file, err := os.Create(fmt.Sprintf("%s.chm", O.inputname))
	if err != nil {
		return err
	}
	defer file.Close()
	_, _ = fmt.Fprint(file, coordline)
	_, _ = fmt.Fprint(file, frozen)
	_, _ = fmt.Fprint(file, "\n\n")
	command := "energy \\\n"
	if Q.Optimize {
		command = "dl-find \\\n"
	}
	_, _ = fmt.Fprint(file, command)
	b := "\\\n             " //break
	sb := "\\\n    "
	link := "0"
	if O.link {
		link = "1"
	}
	mpi:="] \\\n    "
	if O.mpi{
		mpi=fmt.Sprintf("%s mpi_nprocs=$::env(mpi_nprocs) %s mpi_mf=$::env(mpi_mf) %s mpi_omp=$::env(mpi_omp)] %s",b,b,b,sb)
	}
	//I admit the following is horrible. Just take a leap of faith
	arguments := fmt.Sprintf("     theory= %s : [ list basis=%s %s hamiltonian=%s %s accuracy=high %s link=%s %s charge=%d %s jobname=%s %s useghosts=0 %s coords=%s.crd %s %s list_option=full\n\n", O.program, basis, b, method, b, b, link, b, atoms.Charge(), b, O.inputname, b, mpi, O.inputname, sb, optline)
	_, _ = fmt.Fprint(file, arguments)
	_, _ = fmt.Fprint(file, writeline)
	return nonfatal
}

var qcMineDisp = map[string]string{
	"D3":     "d",
	"":       "",
	"nodisp": "",
}

var chemShellMethods = map[string]string{
	"HF":     "hf",
	"hf":     "hf",
	"b3lyp":  "b3-lyp",
	"B3LYP":  "b3-lyp",
	"b3-lyp": "b3-lyp",
	"BLYP":   "blyp",
	"blyp":   "blyp",
	"b-lyp":  "blyp",
	"PBE0":   "pbe0",
	"pbe0":   "pbe0",
	"BP86":   "bp86",
	"bp86":   "bp86",
	"revpbe": "revPBE-d", //qcmine only supports revPBE with -D3 correction
	"revPBE": "revPBE-d",
}

var chemShellBasis = map[string]string{
	"def2-TZVP": "def2-TZVP",
	"tzp":       "def2-TZVP",
	"def2-SVP":  "def2-SVP",
	"svp":       "def2-SVP",
}
