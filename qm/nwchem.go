/*
 * nwchem.go, part of gochem.
 *
 *
 * Copyright 2013 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 "os"
 "strings"
 "fmt"
 "os/exec"
 "github.com/rmera/gochem"
)

//Note that the default methods and basis vary with each program, and even
//for a given program they are NOT considered part of the API, so they can always change.
type NWChemHandle struct {
	defmethod   string
	defbasis    string
	defauxbasis string
	previousMO  string
	restart     bool
	command     string
	inputname   string
	nCPU        int
}

func NewNWChemHandle() *NWChemHandle {
	run := new(NWChemHandle)
	run.SetDefaults()
	return run
}

//NWChemHandle methods

//Sets the number of CPU to be used
func (O *NWChemHandle) SetRestart(r bool) {
	O.restart = r
}

func (O *NWChemHandle) SetName(name string) {
	O.inputname = name
}

func (O *NWChemHandle) SetCommand(name string) {
	O.command = name
}

func (O *NWChemHandle) SetMOName(name string) {
	O.previousMO = name
}

//Sets defaults for NWChem calculation. Default is a single-point at
//TPSS/def2-SVP with RI, and all the available CPU with a max of
//unix.
func (O *NWChemHandle) SetDefaults() {
	O.defmethod = "tpss"
	O.defbasis = "def2-svp"
	O.command =  "nwchem"

}

//BuildInput builds an input for NWChem based int the data in atoms, coords and C.
//returns only error.
func (O *NWChemHandle) BuildInput(coords *chem.VecMatrix, atoms chem.ReadRef, Q *Calc) error {
	//Only error so far
	if atoms == nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
	}
	if Q.Basis == "" {
		fmt.Fprintf(os.Stderr, "no basis set assigned for NWChem calculation, will used the default %s, \n", O.defbasis)
		Q.Basis = O.defbasis
	}
	if Q.Method == "" {
		fmt.Fprintf(os.Stderr, "no method assigned for NWChem calculation, will used the default %s, \n", O.defmethod)
		Q.Method = O.defmethod
		Q.RI = true
	}
	if O.inputname == "" {
		O.inputname = "gochem"
	}
	//The initial guess
	vectors:=fmt.Sprintf("output  %s.movecs",O.inputname)  //The initial guess
	switch Q.Guess{
		case "":
		case "hcore":
			vectors=fmt.Sprintf("input hcore %s",vectors)
		default:
			if !Q.OldMO{
				//If the user gives something in Q.Guess but DOES NOT want an old MO to be used, I assume he/she wants to put whatever 
				//is in Q.Guess directly  in the vector keyword. If you want the default put an empty string in Q.Guess.
				vectors=fmt.Sprintf("%s %s",Q.Guess,vectors)
				break
			}
			//I assume the user gave a basis set name in Q.Guess which I can use to project vectors from a previous run.
			moname:=getOldMO(O.previousMO)
			if moname==""{
				break
			}
			if strings.ToLower(Q.Guess)==strings.ToLower(Q.Basis){
				//Useful if you only change functionals.
				vectors=fmt.Sprintf("input %s %s",moname,vectors)
			}else{
				//This will NOT work if one assigns different basis sets to different atoms.
				vectors=fmt.Sprintf("input project %s %s %s",strings.ToLower(Q.Guess),moname,vectors)
			}
	}
	vectors="vectors "+vectors

	disp,ok := nwchemDisp[Q.Disperssion]
	if !ok{
		disp="vdw 3"
	}

	task := "dft energy"
	driver:=""
	if Q.Optimize == true {
		task = "dft optimize"
		//First an optimization with very loose convergency and the standard trust radius.
		driver=fmt.Sprintf("driver\n maxiter 200\n trust 0.3\n gmax 0.0500\n grms 0.0300\n xmax 0.1800\n xrms 0.1200\n xyz %s_prev\nend\ntask dft optimize",O.inputname)
		//Then the final optimization with a small trust radius and a convergence a bit looser than the nwchem default, which is very tight.
		driver=fmt.Sprintf("%s\ndriver\n maxiter 200\n trust 0.1\n gmax 0.003\n grms 0.0001\n xmax 0.004 \n xrms 0.002\n xyz %s\nend",driver,O.inputname)
	}

	tightness:=""
	switch  Q.SCFTightness {
		case 1:
			tightness="convergence energy 1.000000E-08\n convergence density 5.000000E-09\n convergence gradient 1E-05"
		case 2:
			//NO idea if this will work, or the criteria will be stronger than the criteria for the intergral evaluation
			//and thus the SCF will never converge. Delete when tested.
			tightness="convergence energy 1.000000E-10\n convergence density 5.000000E-11\n convergence gradient 1E-07"
	}

	//Here I dont quite know what to do to help convergency, Ill just slightly extend the iteration tolerance. Sorry about that.
	scfiters:="iterations 30"
	if Q.SCFConvHelp>0{
		scfiters="iterations 60"
	}
	grid,ok:=nwchemGrid[Q.Grid]
	if !ok{
		grid="medium"
	}
	grid=fmt.Sprintf("grid %s",grid)
	var err error

	//Only cartesian constraints supported by now.
	constraints:=""
	if len(Q.CConstraints)>0{
		constraints="constraints\n fix atom"
		for _,v:=range(Q.CConstraints){
			constraints=fmt.Sprintf("%s %d",constraints,v+1) //goChem numbering starts from 0, apparently NWChem starts from 1, hence the v+1
		}
		constraints=constraints+"\nend"
	}

	cosmo := ""
	if Q.Dielectric > 0 {
		if Q.Optimize{
			cosmo=fmt.Sprintf("cosmo\n dielec %4.1f\n do_gasphase True\nend",Q.Dielectric)
		}else{
			cosmo=fmt.Sprintf("cosmo\n dielec %4.1f\n do_gasphase True\nend",Q.Dielectric)
		}
	}
	memory := ""
	if Q.Memory != 0 {
		memory = fmt.Sprintf("memory total %d mb", Q.Memory)
	}
	m:=strings.ToLower(Q.Method)
	method,ok:=nwchemMethods[m]
	if !ok{
		method="xtpss03 ctpss03"
	}
	method=fmt.Sprintf("xc %s",method)

	//////////////////////////////////////////////////////////////
	//Now lets write the thing. Ill process/write the basis later
	//////////////////////////////////////////////////////////////
	file, err := os.Create(fmt.Sprintf("%s.nw", O.inputname))
	if err != nil {
		return err
	}
	defer file.Close()
	start:="start"
	if O.restart{
		start="restart"
	}
	_, err = fmt.Fprintf(file,"%s %s\n",start,O.inputname)
	//after this check its assumed that the file is ok.
	if err != nil {
		return err
	}
	fmt.Fprint(file, "echo\n") //echo input in the output.
	fmt.Fprintf(file, "charge %d\n",atoms.Charge())
	if memory!=""{
		fmt.Fprintf(file,"%s\n", memory) //the memory
	}
	//Now the geometry:
	//If we have cartesian constraints we give the directive noautoz to optimize in cartesian coordinates.
	autoz:=""
	if len(Q.CConstraints)>0{
		autoz="noautoz"
	}
	fmt.Fprintf(file, "geometry units angstroms %s\n",autoz)
	elements:=make([]string,0,5) //I will collect the different elements that are in the molecule using the same loop as the geometry.
	for i := 0; i < atoms.Len(); i++ {
		symbol:=atoms.Atom(i).Symbol
		//In the following if/else I try to set up basis for specific atoms. Not SO sure it works.
		if isInInt(Q.HBAtoms,i){
			symbol=symbol+"1"
		}else if isInInt(Q.LBAtoms,i){
			symbol=symbol+"2"
		}
		fmt.Fprintf(file, " %-2s  %8.3f%8.3f%8.3f \n", symbol, coords.At(i, 0), coords.At(i, 1), coords.At(i, 2))

		if !isInString(elements,symbol){
			elements=append(elements,symbol)
		}
	}
	fmt.Fprintf(file, "end\n")
	//The basis. First the ao basis (required)
	decap:=strings.ToLower //hoping to make the next for loop less ugly

	basis:=make([]string,1,2)
	basis[0]="\"ao basis\""
	fmt.Fprintf(file,"basis \"ao basis\"\n",)
	for _,el:=range(elements){
		if isInString(Q.HBElements,el) || strings.HasSuffix(el,"1"){
			fmt.Fprintf(file," %-2s library %s\n",el,decap(Q.HighBasis))
		}else if isInString(Q.LBElements,el) || strings.HasSuffix(el,"2"){
			fmt.Fprintf(file," %-2s library %s\n",el,decap(Q.LowBasis))
		}else{
			fmt.Fprintf(file, " %-2s library %s\n",el,decap(Q.Basis))
		}
	}
	fmt.Fprintf(file, "end\n")
	//Only Ahlrichs basis are supported for RI. USE AHLRICHS BASIS, PERKELE! :-)
	//The only Ahlrichs J basis in NWchem appear to be equivalent to def2-TZVPP/J (orca nomenclature). I suppose that they are still faster
	//than not using RI if the main basis is SVP. One can also hope that they are good enough if the main basis is QZVPP or something.
	//(about the last point, it appears that in Turbomole, the aux basis also go up to TZVPP).
	//This comment is based on the H, Be and C basis.
	if Q.RI{
		fmt.Fprint(file,"basis \"cd basis\"\n * library \"Ahlrichs Coulomb Fitting\"\nend\n")
	}
	//Now the geometry constraints. I kind of assume they are 
	if constraints!=""{
		fmt.Fprintf(file,"%s\n",constraints)
	}
	if cosmo!=""{
		fmt.Fprintf(file,"%s\n",cosmo)
	}
	//The DFT block
	fmt.Fprint(file,"dft\n")
	fmt.Fprintf(file," %s\n",vectors)
	fmt.Fprintf(file," %s\n",scfiters)
	if tightness!=""{
		fmt.Fprintf(file," %s\n",tightness)
	}
	fmt.Fprintf(file," %s\n",grid)
	fmt.Fprintf(file, " %s\n",method)
	if disp!=""{
		fmt.Fprintf(file," disp %s\n",disp)
	}

	//task part
	fmt.Fprintf(file," mult %d\n",atoms.Multi())
	fmt.Fprint(file,"end\n")
	fmt.Fprintf(file,"%s",driver)
	fmt.Fprintf(file,"task %s\n",task)


	return nil
}







//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//Not waiting for results works
//only for unix-compatible systems, as it uses bash and nohup.
func (O *NWChemHandle) Run(wait bool) (err error) {
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


func getOldMO(prevMO string) string {
	dir, _ := os.Open("./")     //This should always work, hence ignoring the error
	files, _ := dir.Readdir(-1) //Get all the files.
	for _, val := range files {
		if prevMO != "" {
			break
		}
		if val.IsDir() == true {
			continue
		}
		name := val.Name()
		if strings.Contains(".movecs", name) {
			prevMO = name
			break
		}
	}
	if prevMO != "" {
		return prevMO
	} else {
		return ""

	}
}





var nwchemDisp = map[string]string{
	"nodisp": "",
	"D2":     "vdw 2",
	"D3":     "vdw 3",
	"D3ZERO": "vdw 3",
	"D3Zero": "vdw 3",
	"D3zero": "vdw 3",
}

var nwchemGrid = map[int]string{
	1: "xcoarse",
	2: "coarse",
	3: "medium",
	4: "fine",
	5: "xfine",
}


var nwchemMethods = map[string]string{
//	"HF":     "hf",
//	"hf":     "hf",
	"b3lyp":  "b3lyp",
	"b3-lyp": "b3lyp",
//	"PBE":    "pbe",
//	"pbe":    "pbe",
	"pbe0":   "pbe0",
	"revpbe": "revpbe cpbe96",
	"TPSS":   "xtpss03 ctpss03",
	"tpss":   "xtpss03 ctpss03",
	"TPSSh":  "xctpssh",
	"tpssh":  "xctpssh",
	"bp86":   "becke88 perdew86",
	"b-p":    "becke88 perdew86",
	"blyp":   "becke88 lyp",
}






//Reads the latest geometry from an NWChem optimization. Returns the
//geometry or error. Returns the geometry AND error if the geometry read
//is not the product of a correctly ended NWChem calculation. In this case
//the error is "probable problem in calculation"*/
func (O *NWChemHandle) OptimizedGeometry(atoms chem.Ref) (*chem.VecMatrix, error) {
	return nil, fmt.Errorf("not yet implemented")
}

//Gets the energy of a previous NWChem calculation.
//Returns error if problem, and also if the energy returned that is product of an
//abnormally-terminated NWChem calculation. (in this case error is "Probable problem
//in calculation")
func (O *NWChemHandle) Energy() (float64, error) {
	return 0, fmt.Errorf("Not yet implemented")
}


//This checks that an NWChem calculation has terminated normally
//I know this duplicates code, I wrote this one first and then the other one.
func (O *NWChemHandle) nwchemNormalTermination() bool {
	var ini int64 = 0
	var end int64 = 0
	var first bool
	buf := make([]byte, 1)
	f, err := os.Open(fmt.Sprintf("%s.out", O.inputname))
	if err != nil {
		return false
	}
	defer f.Close()
	var i int64 = 1
	for ; ; i++ {
		if _, err := f.Seek(-1*i, 2); err != nil {
			return false
		}
		if _, err := f.Read(buf); err != nil {
			return false
		}
		if buf[0] == byte('\n') && first == false {
			first = true
		} else if buf[0] == byte('\n') && end == 0 {
			end = i
		} else if buf[0] == byte('\n') && ini == 0 {
			ini = i
			break
		}

	}
	f.Seek(-1*ini, 2)
	bufF := make([]byte, ini-end)
	f.Read(bufF)
	if strings.Contains(string(bufF), "**ORCA TERMINATED NORMALLY**") {
		return true
	}
	return false
}
