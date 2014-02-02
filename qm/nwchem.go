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
	"fmt"
	"github.com/rmera/gochem"
	"os"
	"os/exec"
	"strconv"
	"strings"
)

//Note that the default methods and basis vary with each program, and even
//for a given program they are NOT considered part of the API, so they can always change.
type NWChemHandle struct {
	defmethod   string
	defbasis    string
	defauxbasis string
	previousMO  string
	restart     bool
	smartCosmo  bool
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

//Sets the name of a file containing orbitals which will be used as a guess for this calculations
func (O *NWChemHandle) SetMOName(name string) {
	O.previousMO = name
}

//For an optimization, first calculate an SCF with do_gasphase True and use THAT density guess
//for the first optimization step. The optimization is done with do_gasphase False.
//for a SP, smartCosmo simply means do_gasphase False.
//Notice that SmartCosmo is not reallty too smart, for optimizations. In my tests, it doesn't really
//make things better. I keep it for further testing, and may never make it to the master branch.
//My tests indicate that just using do_gasphase False is good enough for optimizations.
func (O *NWChemHandle) SetSmartCosmo(set bool) {
	O.smartCosmo = set
}

//Sets defaults for NWChem calculation. Default is a single-point at
//TPSS/def2-SVP with RI, and all the available CPU with a max of
//unix.
func (O *NWChemHandle) SetDefaults() {
	O.defmethod = "tpss"
	O.defbasis = "def2-svp"
	O.command = "nwchem"
	O.smartCosmo = false

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
	vectors := fmt.Sprintf("output  %s.movecs", O.inputname) //The initial guess
	switch Q.Guess {
	case "":
	case "hcore": //This option is not a great idea, apparently.
		vectors = fmt.Sprintf("input hcore %s", vectors)
	default:
		if !Q.OldMO {
			//If the user gives something in Q.Guess but DOES NOT want an old MO to be used, I assume he/she wants to put whatever
			//is in Q.Guess directly  in the vector keyword. If you want the default put an empty string in Q.Guess.
			vectors = fmt.Sprintf("%s %s", Q.Guess, vectors)
			break
		}
		//I assume the user gave a basis set name in Q.Guess which I can use to project vectors from a previous run.
		moname := getOldMO(O.previousMO)
		if moname == "" {
			break
		}
		if strings.ToLower(Q.Guess) == strings.ToLower(Q.Basis) {
			//Useful if you only change functionals.
			vectors = fmt.Sprintf("input %s %s", moname, vectors)
		} else {
			//This will NOT work if one assigns different basis sets to different atoms.
			vectors = fmt.Sprintf("input project %s %s %s", strings.ToLower(Q.Guess), moname, vectors)
		}
	}
	vectors = "vectors " + vectors

	disp, ok := nwchemDisp[Q.Disperssion]
	if !ok {
		disp = "vdw 3"
	}
	tightness := ""
	switch Q.SCFTightness {
	case 1:
		tightness = "convergence energy 5.000000E-08\n convergence density 5.000000E-09\n convergence gradient 1E-05"
	case 2:
		//NO idea if this will work, or the criteria will be stronger than the criteria for the intergral evaluation
		//and thus the SCF will never converge. Delete when properly tested.
		tightness = "convergence energy 1.000000E-10\n convergence density 5.000000E-11\n convergence gradient 1E-07"
	}

	//For  now, the only convergence help I trust is to run a little HF calculation before and use the orbitals as a guess.
	//It works quite nicely. When the NWChem people get their shit together and fix the bugs with cgmin and RI and cgmin and
	//COSMO, cgmin will be a great option also.
	scfiters := "iterations 60"
	prevscf := ""
	if Q.SCFConvHelp > 0 {
		scfiters = "iterations 200"
		if Q.Guess == "" {
			prevscf = fmt.Sprintf("\nbasis \"3-21g\"\n * library 3-21g\nend\nset \"ao basis\" 3-21g\nscf\n maxiter 200\n vectors output hf.movecs\n %s\nend\ntask scf energy\n\n", strings.ToLower(mopacMultiplicity[atoms.Multi()]))
			vectors = fmt.Sprintf("vectors input project \"3-21g\" hf.movecs output %s.movecs", O.inputname)
		}
	}
	grid, ok := nwchemGrid[Q.Grid]
	if !ok {
		grid = "medium"
	}
	if Q.SCFTightness > 0 { //We need this if we want the SCF to converge.
		grid = "xfine"
	}
	grid = fmt.Sprintf("grid %s", grid)
	var err error

	//Only cartesian constraints supported by now.
	constraints := ""
	if len(Q.CConstraints) > 0 {
		constraints = "constraints\n fix atom"
		for _, v := range Q.CConstraints {
			constraints = fmt.Sprintf("%s %d", constraints, v+1) //goChem numbering starts from 0, apparently NWChem starts from 1, hence the v+1
		}
		constraints = constraints + "\nend"
	}

	cosmo := ""
	if Q.Dielectric > 0 {
		//SmartCosmo in a single-point means that do_gasphase False is used, nothing fancy.
		if Q.Optimize || O.smartCosmo {
			cosmo = fmt.Sprintf("cosmo\n dielec %4.1f\n do_gasphase False\nend", Q.Dielectric)
		} else {
			cosmo = fmt.Sprintf("cosmo\n dielec %4.1f\n do_gasphase True\nend", Q.Dielectric)
		}
	}
	memory := ""
	if Q.Memory != 0 {
		memory = fmt.Sprintf("memory total %d mb", Q.Memory)
	}
	m := strings.ToLower(Q.Method)
	method, ok := nwchemMethods[m]
	if !ok {
		method = "xtpss03 ctpss03"
	}
	method = fmt.Sprintf("xc %s", method)

	task := "dft energy"
	driver := ""
	preopt := ""
	if Q.Optimize == true {
		eprec := "" //The available presition is set to default except if tighter SCF convergene criteria are being used.
		if Q.SCFTightness > 0 {
			eprec = " eprec 1E-7\n"
		}
		if Q.Dielectric > 0 && O.smartCosmo {
			//If COSMO is used, and O.SmartCosmo is enabled, we start the optimization with a rather loose SCF (the default).
			//and use that denisty as a starting point for the next calculation. The idea is to
			//avoid the gas phase calculation in COSMO.
			//This procedure doesn't seem to help at all, and just using do_gasphase False appears to be good enough in my tests.
			preopt = fmt.Sprintf("cosmo\n dielec %4.1f\n do_gasphase True\nend\n", Q.Dielectric)
			preopt = fmt.Sprintf("%sdft\n iterations 100\n %s\n %s\n print low\nend\ntask dft energy\n", preopt, vectors, method)
			vectors = fmt.Sprintf("vectors input %s.movecs output  %s.movecs", O.inputname, O.inputname) //We must modify the initial guess so we use the vectors we have just generated
		}
		//The NWCHem optimizer is horrible. To try to get something out of it we use this 3-step optimization scheme where we try to compensate for the lack of
		//variable trust radius in nwchem.
		task = "dft optimize"
		//First an optimization with very loose convergency and the standard trust radius.
		driver = fmt.Sprintf("driver\n maxiter 200\n%s trust 0.3\n gmax 0.0500\n grms 0.0300\n xmax 0.1800\n xrms 0.1200\n xyz %s_prev\nend\ntask dft optimize", eprec, O.inputname)
		//Then a second optimization with a looser convergency and a 0.1 trust radius
		driver = fmt.Sprintf("%s\ndriver\n maxiter 200\n%s trust 0.1\n gmax 0.009\n grms 0.001\n xmax 0.04 \n xrms 0.02\n xyz %s_prev2\nend\ntask dft optimize", driver, eprec, O.inputname)
		//Then the final optimization with a small trust radius and the NWChem default convergence criteria.
		driver = fmt.Sprintf("%s\ndriver\n maxiter 200\n%s trust 0.05\n xyz %s\nend\n", driver, eprec, O.inputname)
		//Old criteria (ORCA): gmax 0.003\n grms 0.0001\n xmax 0.004 \n xrms 0.002\n
	}

	//////////////////////////////////////////////////////////////
	//Now lets write the thing. Ill process/write the basis later
	//////////////////////////////////////////////////////////////
	file, err := os.Create(fmt.Sprintf("%s.nw", O.inputname))
	if err != nil {
		return err
	}
	defer file.Close()
	start := "start"
	if O.restart {
		start = "restart"
	}
	_, err = fmt.Fprintf(file, "%s %s\n", start, O.inputname)
	//after this check its assumed that the file is ok.
	if err != nil {
		return err
	}
	fmt.Fprint(file, "echo\n") //echo input in the output.
	fmt.Fprintf(file, "charge %d\n", atoms.Charge())
	if memory != "" {
		fmt.Fprintf(file, "%s\n", memory) //the memory
	}
	//Now the geometry:
	//If we have cartesian constraints we give the directive noautoz to optimize in cartesian coordinates.
	autoz := ""
	if len(Q.CConstraints) > 0 {
		autoz = "noautoz"
	}
	fmt.Fprintf(file, "geometry units angstroms noautosym %s\n", autoz)
	elements := make([]string, 0, 5) //I will collect the different elements that are in the molecule using the same loop as the geometry.
	for i := 0; i < atoms.Len(); i++ {
		symbol := atoms.Atom(i).Symbol
		//In the following if/else I try to set up basis for specific atoms. Not SO sure it works.
		if isInInt(Q.HBAtoms, i) {
			symbol = symbol + "1"
		} else if isInInt(Q.LBAtoms, i) {
			symbol = symbol + "2"
		}
		fmt.Fprintf(file, " %-2s  %8.3f%8.3f%8.3f \n", symbol, coords.At(i, 0), coords.At(i, 1), coords.At(i, 2))

		if !isInString(elements, symbol) {
			elements = append(elements, symbol)
		}
	}
	fmt.Fprintf(file, "end\n")
	fmt.Fprintf(file, prevscf) //The preeliminar SCF if exists.
	//The basis. First the ao basis (required)
	decap := strings.ToLower //hoping to make the next for loop less ugly
	basis := make([]string, 1, 2)
	basis[0] = "\"ao basis\""
	fmt.Fprintf(file, "basis \"large\" spherical\n") //According to the manual this fails with COSMO. The calculations dont crash. Need to compare energies and geometries with Turbomole in order to be sure.
	for _, el := range elements {
		if isInString(Q.HBElements, el) || strings.HasSuffix(el, "1") {
			fmt.Fprintf(file, " %-2s library %s\n", el, decap(Q.HighBasis))
		} else if isInString(Q.LBElements, el) || strings.HasSuffix(el, "2") {
			fmt.Fprintf(file, " %-2s library %s\n", el, decap(Q.LowBasis))
		} else {
			fmt.Fprintf(file, " %-2s library %s\n", el, decap(Q.Basis))
		}
	}
	fmt.Fprintf(file, "end\n")
	fmt.Fprintf(file, "set \"ao basis\" large\n")
	//Only Ahlrichs basis are supported for RI. USE AHLRICHS BASIS, PERKELE! :-)
	//The only Ahlrichs J basis in NWchem appear to be equivalent to def2-TZVPP/J (orca nomenclature). I suppose that they are still faster
	//than not using RI if the main basis is SVP. One can also hope that they are good enough if the main basis is QZVPP or something.
	//(about the last point, it appears that in Turbomole, the aux basis also go up to TZVPP).
	//This comment is based on the H, Be and C basis.
	if Q.RI {
		fmt.Fprint(file, "basis \"cd basis\"\n * library \"Ahlrichs Coulomb Fitting\"\nend\n")
	}
	//Now the geometry constraints. I kind of assume they are
	if constraints != "" {
		fmt.Fprintf(file, "%s\n", constraints)
	}
	fmt.Fprintf(file, preopt)
	if cosmo != "" {
		fmt.Fprintf(file, "%s\n", cosmo)
	}
	//The DFT block
	fmt.Fprint(file, "dft\n")
	fmt.Fprintf(file, " %s\n", vectors)
	fmt.Fprintf(file, " %s\n", scfiters)
	if tightness != "" {
		fmt.Fprintf(file, " %s\n", tightness)
	}
	fmt.Fprintf(file, " %s\n", grid)
	fmt.Fprintf(file, " %s\n", method)
	if disp != "" {
		fmt.Fprintf(file, " disp %s\n", disp)
	}
	if Q.Optimize {
		fmt.Fprintf(file, " print convergence\n")
	}
	//task part
	fmt.Fprintf(file, " mult %d\n", atoms.Multi())
	fmt.Fprint(file, "end\n")
	fmt.Fprintf(file, "%s", driver)
	fmt.Fprintf(file, "task %s\n", task)

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
	"b3lyp":   "b3lyp",
	"b3-lyp":  "b3lyp",
	"pbe0":    "pbe0",
	"mpw1b95": "mpw1b95",
	"revpbe":  "revpbe cpbe96",
	"TPSS":    "xtpss03 ctpss03",
	"tpss":    "xtpss03 ctpss03",
	"TPSSh":   "xctpssh",
	"tpssh":   "xctpssh",
	"bp86":    "becke88 perdew86",
	"b-p":     "becke88 perdew86",
	"blyp":    "becke88 lyp",
}

//Reads the latest geometry from an NWChem optimization. Returns the
//geometry or error. Returns the geometry AND error if the geometry read
//is not the product of a correctly ended NWChem calculation. In this case
//the error is "probable problem in calculation".
func (O *NWChemHandle) OptimizedGeometry(atoms chem.Ref) (*chem.VecMatrix, error) {
	var err2 error
	lastnumber := 0
	lastname := ""
	if !O.nwchemNormalTermination() {
		err2=fmt.Errorf("Probable problem in calculation")
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
	if err!=nil{
		return nil, err
	}
	return mol.Coords[0], err2
}

//Gets the energy of a previous NWChem calculation.
//Returns error if problem, and also if the energy returned that is product of an
//abnormally-terminated NWChem calculation. (in this case error is "Probable problem
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

//This checks that an NWChem calculation has terminated normally
//I know this duplicates code, I wrote this one first and then the other one.
func (O *NWChemHandle) nwchemNormalTermination() bool {
	ret:=false
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
			ret=true
			break
		}
	}
	return ret
}

