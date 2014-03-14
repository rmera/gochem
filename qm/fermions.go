/*
 * fermions.go, part of gochem.
 *
 *
 * Copyright 2014 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
	"github.com/rmera/gochem"
	"os"
//	"os/exec"
//	"strconv"
	"strings"
)

//Note that the default methods and basis vary with each program, and even
//for a given program they are NOT considered part of the API, so they can always change.
type FermionsHandle struct {
	defmethod   string
	defbasis    string
	defauxbasis string
//	previousMO  string
//	restart     bool
//	smartCosmo  bool
	command     string
	inputname   string
	nCPU        int
	gpu         string
}

func NewFermionsHandle() *FermionsHandle {
	run := new(FermionsHandle)
	run.SetDefaults()
	return run
}

//FermionsHandle methods


//Sets GPU usage. Alternatives are "cuda" or "opencl" (alias "ocl"). Anything else is ignored. GPU is off by default.
func (O *FermionsHandle) SetGPU(rawname string){
	name:=strings.ToLower(rawname)
	if name=="cuda"{
		O.gpu="USE_CUDA YES\nCUDA_LINALG_MIN 500"
	}else if name=="opencl" || name=="ocl"{
		O.gpu="USE_OCL YES\nOCL_LINALG_MIN 500" //OpenCL is only for SCF energies!!!!!
	}
}

func (O *FermionsHandle) SetName(name string) {
	O.inputname = name
}

func (O *FermionsHandle) SetCommand(name string) {
	O.command = name
}


//Sets defaults for Fermions++ calculation. Default is a single-point at
//revPBE/def2-SVP
func (O *FermionsHandle) SetDefaults() {
	O.defmethod = "revpbe"
	O.defbasis = "def2-svp"
	O.command = ""
	O.inputname = "gochem"
	O.gpu=""

}

//BuildInput builds an input for Fermions++ based int the data in atoms, coords and C.
//returns only error.
func (O *FermionsHandle) BuildInput(coords *chem.VecMatrix, atoms chem.ReadRef, Q *Calc) error {
	//Only error so far


	if atoms == nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
	}
	if Q.Basis == "" {
		fmt.Fprintf(os.Stderr, "no basis set assigned for Fermions++ calculation, will used the default %s, \n", O.defbasis)
		Q.Basis = O.defbasis
	}
	if Q.Method == "" {
		fmt.Fprintf(os.Stderr, "no method assigned for Fermions++ calculation, will used the default %s, \n", O.defmethod)
		Q.Method = O.defmethod
	}

	disp,ok:=fermionsDisp[strings.ToLower(Q.Disperssion)]
	if !ok{
		disp = "disp_corr  D3"
	}

	grid, ok :=fermionsGrid[Q.Grid]
	if !ok {
		grid = "M3"
	}
	grid = fmt.Sprintf("GRID_RAD_TYPE %s", grid)
	var err error

	m := strings.ToLower(Q.Method)
	method, ok := fermionsMethods[m]
	if !ok {
		method = "EXC XC_GGA_X_PBE_R\n ECORR XC_GGA_C_PBE"
	}
	task := "SinglePoint"
	dloptions:=""
	if Q.Optimize == true {
		task="DLF_OPTIMIZE"
		dloptions=fmt.Sprintf("*start::dlfind\n JOB std\n method l-bfgs\n trust_radius energy\n dcd %s.dcd\n maxcycle 300\n maxene 200\n coord_type cartesian\n*end\n",O.inputname)
		//Only cartesian constraints supported by now.
		if len(Q.CConstraints) > 0 {
			dloptions=fmt.Sprintf("%s\n*start::dlconstraints\n",dloptions)
			for _,v:=range(Q.CConstraints){
				dloptions=fmt.Sprintf("%s cart %d\n",dloptions,v+1) //fortran numbering, starts with 1.
			}
			dloptions=fmt.Sprintf("%s*end\n",dloptions)
		}
	}
	cosmo := ""
	if Q.Dielectric > 0 {
		cosmo = fmt.Sprintf("*start::solvate\n pcm_model cpcm\n epsilon %f\n cavity_model bondi\n*end\n", Q.Dielectric)
	}

	//////////////////////////////////////////////////////////////
	//Now lets write the thing.
	//////////////////////////////////////////////////////////////
	file, err := os.Create(fmt.Sprintf("%s.in", O.inputname))
	if err != nil {
		return err
	}
	defer file.Close()
	//Start with the geometry part (coords, charge and multiplicity)
	fmt.Fprintf(file,"*start::geo\n")
	fmt.Fprintf(file, "%d %d\n", atoms.Charge(), atoms.Multi())
	for i := 0; i < atoms.Len(); i++ {
		fmt.Fprintf(file, "%-2s  %8.3f%8.3f%8.3f\n", atoms.Atom(i).Symbol, coords.At(i, 0), coords.At(i, 1), coords.At(i, 2))
	}
	fmt.Fprintf(file,"*end\n\n")
	fmt.Fprintf(file,"*start::sys\n")
	fmt.Fprintf(file," TODO %s\n",task)
	fmt.Fprintf(file," BASIS %s\n PC pure\n",strings.ToLower(Q.Basis))
	fmt.Fprintf(file," %s\n",method)
	fmt.Fprintf(file," %s\n",grid)
	fmt.Fprintf(file," %s\n",disp)
	fmt.Fprintf(file," %s\n",O.gpu)
	if !Q.Optimize{
		fmt.Fprintf(file," INFO 2\n")
	}
	fmt.Fprintf(file,"*end\n\n")
	fmt.Fprintf(file,"%s\n",cosmo)
	fmt.Fprintf(file,"%s\n",dloptions)
	return nil
}

/*********
//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//Not waiting for results works
//only for unix-compatible systems, as it uses bash and nohup.
func (O *FermionsHandle) Run(wait bool) (err error) {
	if wait == true {
		out, err := os.Create(fmt.Sprintf("%s.out", O.inputname))
		if err != nil {
			return err
		}
		defer out.Close()
		command := exec.Command(O.command, fmt.Sprintf("%s.nw", O.inputname))
		command.Stdout = out
		err = command.Run()

	} else {
		//This will not work in windows.
		command := exec.Command("sh", "-c", "nohup "+O.command+fmt.Sprintf(" %s.nw > %s.out &", O.inputname, O.inputname))
		err = command.Start()
	}
	return err
}
*******/

var fermionsDisp = map[string]string{
	"":       "",
	"nodisp": "",
	"d2":     "disp_corr D2",
	"d3bj":   "disp_corr BJ",
	"bj":     "disp_corr BJ",
	"d3":     "disp_corr D3",
}


//M5 etc are not supported
var fermionsGrid = map[int]string{
	1: "M3",
	2: "M3",
	3: "M3",
	4: "M4",
	5: "M4",
}

//All meta-gga functionals here are buggy in the libxc library, and thus not reliable. I leave them hoping this will get fixed eventually.
var fermionsMethods = map[string]string{
	"b3lyp":   "EXC XC_HYB_GGA_XC_B3LYP\n ECORR XC_HYB_GGA_XC_B3LYP",
	"b3-lyp":   "EXC XC_HYB_GGA_XC_B3LYP\n ECORR XC_HYB_GGA_XC_B3LYP",
	"pbe":      "EXC XC_GGA_X_PBE\n ECORR XC_GGA_C_PBE",
	"revpbe":   "EXC XC_GGA_X_PBE_R\n ECORR XC_GGA_C_PBE",
	"pbe0":    "EXC XC_HYB_GGA_XC_PBEH\n ECORR XC_HYB_GGA_XC_PBEH",
	"mpw1b95": "EXC XC_HYB_MGGA_XC_MPW1B95\n ECORR XC_HYB_MGGA_XC_MPW1B95", //probably also buggy
	"tpss":    "EXC XC_MGGA_X_TPSS\n ECORR XC_MGGA_C_TPSS",  //Buggy in current libxc, avoid.
	"bp86":    "EXC XC_GGA_X_B88\n ECORR XC_GGA_C_P86",
	"b-p":     "EXC XC_GGA_X_B88\n ECORR XC_GGA_C_P86",
	"blyp":    "EXC XC_GGA_X_B88\n ECORR XC_GGA_C_LYP",
	"b-lyp":    "EXC XC_GGA_X_B88\n ECORR XC_GGA_C_LYP",
}

/****************************************************
//Reads the latest geometry from an Fermions++ optimization. Returns the
//geometry or error. Returns the geometry AND error if the geometry read
//is not the product of a correctly ended Fermions++ calculation. In this case
//the error is "probable problem in calculation".
func (O *FermionsHandle) OptimizedGeometry(atoms chem.Ref) (*chem.VecMatrix, error) {
	var err2 error
	lastnumber := 0
	lastname := ""
	if !O.nwchemNormalTermination() {
		err2 = fmt.Errorf("Probable problem in calculation")
	}
	dir, err := os.Open("./")
	if err != nil {
		return nil, err
	}
	files, err := dir.Readdirnames(-1)
	if err != nil {
		return nil, err
	}
	//This is a crappy sort/filter, but really, it will never be the bottleneck.
	//We go over the dir content and look for xyz files with the name of the input and without
	// the susbstring _prev. Among these, we choose the file with the largest number in the filename
	//(which will be the latest geometry written) and return that geometry.
	for _, v := range files {
		//quite ugly, feel free to fix.
		if !(strings.Contains(v, O.inputname) && strings.Contains(v, ".xyz") && !strings.Contains(v, "_prev")) {
			continue
		}
		numbers := strings.Split(strings.Replace(v, ".xyz", "", 1), "-")
		ndx, err := strconv.Atoi(numbers[len(numbers)-1])
		if err != nil {
			continue
		}
		if ndx >= lastnumber {
			lastnumber = ndx
			lastname = v
		}
	}
	if lastname == "" {
		return nil, fmt.Errorf("Geometry not found")
	}
	mol, err := chem.XYZRead(lastname)
	if err != nil {
		return nil, err
	}
	return mol.Coords[0], err2
}

//Gets the energy of a previous Fermions++ calculation.
//Returns error if problem, and also if the energy returned that is product of an
//abnormally-terminated Fermions++ calculation. (in this case error is "Probable problem
//in calculation")
func (O *FermionsHandle) Energy() (float64, error) {
	err := fmt.Errorf("Probable problem in calculation")
	f, err1 := os.Open(fmt.Sprintf("%s.out", O.inputname))
	if err1 != nil {
		return 0, err1
	}
	defer f.Close()
	f.Seek(-1, 2) //We start at the end of the file
	energy := 0.0
	var found bool
	for i := 0; ; i++ {
		line, err1 := getTailLine(f)
		if err1 != nil {
			return 0.0, err1
		}
		if strings.Contains(line, "CITATION") {
			err = nil
		}
		if strings.Contains(line, "Total DFT energy") {
			splitted := strings.Fields(line)
			energy, err1 = strconv.ParseFloat(splitted[len(splitted)-1], 64)
			if err1 != nil {
				return 0.0, err1
			}
			found = true
			break
		}
	}
	if !found {
		return 0.0, fmt.Errorf("Output does not contain energy")
	}
	return energy * chem.H2Kcal, err
}

//This checks that an Fermions++ calculation has terminated normally
//I know this duplicates code, I wrote this one first and then the other one.
func (O *FermionsHandle) nwchemNormalTermination() bool {
	ret := false
	f, err1 := os.Open(fmt.Sprintf("%s.out", O.inputname))
	if err1 != nil {
		return false
	}
	defer f.Close()
	f.Seek(-1, 2) //We start at the end of the file
	for i := 0; ; i++ {
		line, err1 := getTailLine(f)
		if err1 != nil {
			return false
		}
		if strings.Contains(line, "CITATION") {
			ret = true
			break
		}
	}
	return ret
}

**********/
