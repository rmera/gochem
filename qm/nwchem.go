/*
 * qm.go, part of gochem.
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

import "os"
import "strings"
import "fmt"
//import "runtime"
import "os/exec"
import "strconv"
import "github.com/rmera/gochem"

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
func (O *NWChemHandle) SetnCPU(cpu int) {
	O.nCPU = cpu
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

/*Sets defaults for ORCA calculation. Default is a single-point at
revPBE/def2-SVP with RI, and all the available CPU with a max of
8. The ORCA command is set to $ORCA_PATH/orca, at least in
unix.*/
func (O *NWChemHandle) SetDefaults() {
	O.defmethod = "tpss"
	O.defbasis = "def2-svp"
//	O.defauxbasis = "def2-SVP/J"
	O.command =  "nw"   //os.ExpandEnv("${ORCA_PATH}/orca")

//	cpu := runtime.NumCPU()
//	if cpu > 8 {
//		O.nCPU = 8
//	}

}

//BuildInput builds an input for ORCA based int the data in atoms, coords and C.
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
			if Q.Guess==Q.Basis{
				//Useful if you only change functionals.
				vectors=fmt.Sprintf("input %s %s",moname,vectors)
			}else{
				//This will NOT work if one assigns different basis sets to different atoms.
				vectors=fmt.Sprintf("input project %s %s %s",Q.Guess,moname,vectors)
			}
		
	}

	disp,ok := nwchemDisp[Q.Disperssion]
	if !ok{
		disp="vdw 3"
	}
	task := "dft energy"
	if Q.Optimize == true {
		task = "dft optimize"
	}

	tightness:=""
	switch  Q.SCFTightness {
		case 1:
			tightness:="convergence energy 1.000000E-08\nconvergence density 5.000000E-09\nconvergence gradient 1E-05"
		case 2:
			//NO idea if this will work, or the criteria will be stronger than the criteria for the intergral evaluation
			//and thus the SCF will never converge. Delete when tested.
			tightness="convergence energy 1.000000E-10\nconvergence density 5.000000E-11\nconvergence gradient 1E-07"
	}

	//Here I dont quite know what to do to help convergency, Ill just slightly extend the iteration tolerance. Sorry about that.
	scfiters:="iterations 30"
	if Q.SCFConvHelp>0{
		scfiters="iterations 50"
	}
	grid,ok:=nwchemGrid[Q.Grid]
	if !ok{
		grid="m"
	}
	var err error

	//Only cartesian constraints supported by now.
	constraints:=""
	if len(Q.CConstraints)>0{
		constraints="constraints\n fix atom"
		for _,v:=range(Q.CConstraints){
			constraints=fmt.Sprintf("%s %i",constraints,v+1) //goChem numbering starts from 0, apparently NWChem starts from 1, hence the v+1
		}
		constraints=constraints+"\nend"
	}

	cosmo := ""
	if Q.Dielectric > 0 {
		cosmo=fmt.Sprintf("cosmo\n dielec %s\nend")
	}
	mem := ""
	if Q.Memory != 0 {
		mem = fmt.Sprintf("memory %d mb", Q.Memory)
	}



	//////////////////////////////////////////////////////////////
	//Now lets write the thing. Ill process/write the basis later
	//////////////////////////////////////////////////////////////
	if O.inputname == "" {
		O.inputname = "gochem"
	}
	file, err := os.Create(fmt.Sprintf("%s.nw", O.inputname))
	if err != nil {
		return err
	}
	defer file.Close()
	start:="start"
	if O.restart{
		start="restart"
	}
	_, err = fmt.Fprint(file,"%s %s\n",start,O.inputname)
	//after this check its assumed that the file is ok.
	if err != nil {
		return err
	}
	fmt.Fprint(file, "echo\n") //echo input in the output.
	fmt.Fprintf(file, "charge %d",atoms.Charge())
	fmt.Fprintf(file,"%s\n", mem) //the memory

	//Now the geometry:
	fmt.Fprint(file, "geometry units angstroms\n")
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
	for _,el:=range(elements){
		if isInString(Q.HBElements,el) || strings.HasSuffix(el,"1"){
			fmt.Fprintf(file,"%s library %s",el,HighBasis)
		}else if isInString(Q.LBElements,el) || strings.HasSuffix(el,"2"){
			fmt.Fprintf(file,"%s library %s",el,LowBasis)
		}else{
			fmt.Fprintf(fie, "%s library %s",el,Basis)
		}
	}
	fmt.Fprint(file, mem)
	fmt.Fprint(file, constraints)
	fmt.Fprint(file, iconstraints)
	fmt.Fprint(file, ElementBasis)
	fmt.Fprint(file, cosmo)
	fmt.Fprint(file, "\n")
	//Now the type of coords, charge and multiplicity
	fmt.Fprintf(file, "* xyz %d %d\n", atoms.Charge(), atoms.Multi())
	//now the coordinates
	//	fmt.Println(atoms.Len(), coords.Rows()) ///////////////

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
		command := exec.Command(O.command, fmt.Sprintf("%s.inp", O.inputname))
		command.Stdout = out
		err = command.Run()

	} else {
		command := exec.Command("sh", "-c", "nohup "+O.command+fmt.Sprintf(" %s.inp > %s.out &", O.inputname, O.inputname))
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











//buildIConstraints transforms the list of cartesian constrains in the QMCalc structre
//into a string with ORCA-formatted internal constraints.
func (O *NWChemHandle) buildIConstraints(C []*IConstraint) (string, error) {
	if C == nil {
		return "\n", nil //no constraints
	}
	constraints := make([]string, len(C)+3)
	constraints[0] = "%geom Constraints\n"
	for key, val := range C {

		if iConstraintOrder[val.Class] != len(val.CAtoms) {
			return "", fmt.Errorf("Internal constraint ill-formated")
		}

		var temp string
		if val.Class == 'B' {
			temp = fmt.Sprintf("         {B %d %d %2.3f C}\n", val.CAtoms[0], val.CAtoms[1], val.Val)
		} else if val.Class == 'A' {
			temp = fmt.Sprintf("         {A %d %d %d %2.3f C}\n", val.CAtoms[0], val.CAtoms[1], val.CAtoms[2], val.Val)
		} else if val.Class == 'D' {
			temp = fmt.Sprintf("         {D %d %d %d %d %2.3f C}\n", val.CAtoms[0], val.CAtoms[1], val.CAtoms[2], val.CAtoms[3], val.Val)
		}
		constraints[key+1] = temp
	}
	last := len(constraints) - 1
	constraints[last-1] = "         end\n"
	constraints[last] = " end\n"
	final := strings.Join(constraints, "")
	return final, nil
}

var iConstraintOrder = map[byte]int{
	'B': 2,
	'A': 3,
	'D': 4,
}

//buildCConstraints transforms the list of cartesian constrains in the QMCalc structre
//into a string with ORCA-formatted cartesian constraints
func (O *NWChemHandle) buildCConstraints(C []int) string {
	if C == nil {
		return "\n" //no constraints
	}
	constraints := make([]string, len(C)+3)
	constraints[0] = "%geom Constraints\n"
	for key, val := range C {
		constraints[key+1] = fmt.Sprintf("         {C %d C}\n", val)
	}
	last := len(constraints) - 1
	constraints[last-1] = "         end\n"
	constraints[last] = " end\n"
	final := strings.Join(constraints, "")
	return final
}


var orcaSCFTight = map[int]string{
	0: "",
	1: "TightSCF",
	2: "VeryTightSCF",
}

var orcaSCFConv = map[int]string{
	0: "",
	1: "SlowConv",
	2: "VerySlowConv",
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
	"B3LYP":  "b3lyp",
	"b3-lyp": "b3lyp",
//	"PBE":    "pbe",
//	"pbe":    "pbe",
	"pbe0":   "pbe0",
	"PBE0":   "pbe0",
	"TPSS":   "xtpss03 ctpss03",
	"tpss":   "xtpss03 ctpss03",
	"TPSSh":  "xctpssh",
	"tpssh":  "xctpssh",
	"BP86":   "becke88 perdew 86",
	"b-p":    "becke88 perdew 86",
	"blyp":   "becke88 lyp",
}






/*Reads the latest geometry from an ORCA optimization. Returns the
  geometry or error. Returns the geometry AND error if the geometry read
  is not the product of a correctly ended ORCA calculation. In this case
  the error is "probable problem in calculation"*/
func (O *NWChemHandle) OptimizedGeometry(atoms chem.Ref) (*chem.VecMatrix, error) {
	var err error
	geofile := fmt.Sprintf("%s.xyz", O.inputname)
	//Here any error of orcaNormal... or false means the same, so the error can be ignored.
	if trust := O.orcaNormalTermination(); !trust {
		err = fmt.Errorf("Probable problem in calculation")
	}
	//This might not be super efficient but oh well.
	mol, err1 := chem.XYZRead(geofile)
	if err1 != nil {
		return nil, err1
	}
	return mol.Coords[0], err //returns the coords, the error indicates whether the structure is trusty (normal calculation) or not
}

//Gets the energy of a previous Orca calculations.
//Returns error if problem, and also if the energy returned that is product of an
//abnormally-terminated ORCA calculation. (in this case error is "Probable problem
//in calculation")
func (O *NWChemHandle) Energy() (float64, error) {
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
		if strings.Contains(line, "**ORCA TERMINATED NORMALLY**") {
			err = nil
		}
		if strings.Contains(line, "FINAL SINGLE POINT ENERGY") {
			splitted := strings.Fields(line)
			energy, err1 = strconv.ParseFloat(splitted[4], 64)
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

//Gets previous line of the file f
func getTailLine(f *os.File) (line string, err error) {
	var i int64 = 1
	buf := make([]byte, 1)
	var ini int64 = -1
	for ; ; i++ {
		//move the pointer back one byte per cycle
		if _, err := f.Seek(-2, 1); err != nil {
			return "", err
		}
		if _, err := f.Read(buf); err != nil {
			return "", err
		}
		if buf[0] == byte('\n') && ini == -1 {
			ini = i
			break
		}
	}
	bufF := make([]byte, ini)
	f.Read(bufF)

	if _, err := f.Seek(int64(-1*(len(bufF))), 1); err != nil { //making up for the read
		return "", err
	}
	return string(bufF), nil
}

//This checks that an ORCA calculation has terminated normally
//I know this duplicates code, I wrote this one first and then the other one.
func (O *NWChemHandle) orcaNormalTermination() bool {
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
