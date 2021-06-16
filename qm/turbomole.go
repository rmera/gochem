/*
 * turbomole.go, part of gochem.
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


//The TM handler implementation differs from the rest in that it uses several TM programs
// (define, x2t, t2x, cosmoprep) in order to prepare the input and retrieve results.
//Because of this, the programs using this handler will not work if TM is not installed.
//The handler has been made to work with TM7.

package qm

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)


//TMHandle is the representation of a Turbomole (TM) calculation
//This imlpementation supports only singlets and doublets.
type TMHandle struct {
	defmethod   string
	defbasis    string
	defauxbasis string
	previousMO  string
	command     string
	inputname   string
	gimic       bool
	marij       bool
	dryrun      bool
}

//Creates and initializes a new instance of TMRuner, with values set
//to its defaults.
func NewTMHandle() *TMHandle {
	run := new(TMHandle)
	run.SetDefaults()
	return run
}

const noCosmoPrep = "goChem/QM: Unable to run cosmoprep"

//TMHandle methods

//SetName sets the name of the subdirectory, in the current directory
//where the calculation will be ran
func (O *TMHandle) SetName(name string) {
	O.inputname = name

}

//SetMARIJ sets the multipole acceleration
func (O *TMHandle) SetMARIJ(state bool) {
	O.marij = state
}

//SetDryRun sets the flag to see this is a dry run or
//if define will actually be run.
func (O *TMHandle) SetDryRun(dry bool) {
	O.dryrun = dry
}

//SetCommand doesn't do anything, and it is here only for compatibility.
//In TM the command is set according to the method. goChem assumes a normal TM installation.
func (O *TMHandle) SetCommand(name string) {
	//Does nothing again
}

//SetDefaults sets default values for TMHandle. default is an optimization at
//  TPSS-D3 / def2-SVP
//Defaults are not part of the API, they might change as new methods appear.
func (O *TMHandle) SetDefaults() {
	O.defmethod = "tpss"
	O.defbasis = "def2-SVP"
	O.defauxbasis = "def2-SVP"
	O.command = "ridft"
	O.marij = false  //Apparently marij can cause convergence problems
	O.dryrun = false //define IS run by default.
	O.inputname = "gochemturbo"
}

//addMARIJ adds the multipole acceleration if certain conditions are fullfilled:
//O.marij must be true
//The RI approximation must be in use
//The system must have more than 20 atoms
//The basis set cannot be very large (i.e. it can NOT be quadruple-zeta, tzvpp, or basis with diffuse functions)
func (O *TMHandle) addMARIJ(defstring string, atoms chem.AtomMultiCharger, Q *Calc) string {
	if !O.marij {
		return defstring
	}
	if strings.Contains(strings.ToLower(Q.Basis), "def2") && strings.HasSuffix(Q.Basis, "d") { //Rappoport basis
		return defstring
	}
	if strings.Contains(Q.Basis, "cc") && strings.Contains(Q.Basis, "aug") { //correlation consistent with diffuse funcs.
		return defstring
	}
	if strings.Contains(strings.ToLower(Q.Basis), "qz") { //both cc-pVQZ and def2-QZVP and QZVPP
		return defstring
	}

	if strings.Contains(strings.ToLower(Q.Basis), "def2-tzvpp") { //This is less clear but just in case I won't add MARIJ for tzvpp
		return defstring
	}
	//I Have no idea what the MARIJ string was supposed to be. For now the method doesn't work. Remove this if you fix the code
	//below, and uncoment it.
	//	if Q.RI && atoms.Len() >= 20 {
	//		defstring = fmt.Sprintf("%s%s\n\n")
	//	}
	return defstring
}

//Adds all the strings in toapend to the control file, just before the $symmetry keyword
func (O *TMHandle) addToControl(toappend []string, Q *Calc) error {
	f, err := os.Open("control")
	if err != nil {
		return Error{ErrCantInput, Turbomole, O.inputname, "", []string{"os.Open", "addtoControl"}, true}
	}
	lines := make([]string, 0, 200) //200 is just a guess for the number of lines in the control file
	c := bufio.NewReader(f)
	for err == nil {
		var line string
		line, err = c.ReadString('\n')
		lines = append(lines, line)
	}
	f.Close() //I cant defer it because I need it closed now.
	out, err := os.Create("control")
	if err != nil {
		return Error{ErrCantInput, Turbomole, O.inputname, "", []string{"os.Create", "addtoControl"}, true}
	}
	defer out.Close()
	var k string
	for _, i := range lines {
		k = i //May not be too efficient
		if strings.Contains(i, "$symmetry") {
			for _, j := range toappend {
				if _, err := fmt.Fprintf(out, j+"\n"); err != nil {
					return Error{ErrCantInput, Turbomole, O.inputname, "", []string{"fmt.Fprintf", "addtoControl"}, true}
				}
			}
		}
		if Q.SCFConvHelp >= 1 {
			if strings.Contains(k, "$scfiterlimit") {
				k = "$scfiterlimit   100\n"
			}
			if strings.Contains(k, "$scfdamp") {
				k = "$scfdamp start=10 step=0.005 min=0.5\n"
			}
		}
		if _, err := fmt.Fprintf(out, k); err != nil {
			return Error{ErrCantInput, Turbomole, O.inputname, "", []string{"fmt.Fprintf", "addtoControl"}, true}

		}
	}
	return nil
}

func (O *TMHandle) addCosmo(epsilon float64) error {
	//The ammount of newlines is wrong, must fix
	cosmostring := "" //a few newlines before the epsilon
	if epsilon == 0 {
		return nil
	}
	cosmostring = fmt.Sprintf("%s%3.1f\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nr all b\n*\n\n\n\n\n\n", cosmostring, epsilon)
	def := exec.Command("cosmoprep")
	pipe, err := def.StdinPipe()
	if err != nil {
		return Error{noCosmoPrep, Turbomole, O.inputname, err.Error(), []string{"exec.StdinPipe", "addCosmo"}, true}
	}
	defer pipe.Close()
	pipe.Write([]byte(cosmostring))
	if err := def.Run(); err != nil {
		return Error{noCosmoPrep, Turbomole, O.inputname, err.Error(), []string{"exec.Run", "addCosmo"}, true}

	}
	return nil

}

func (O *TMHandle) addBasis(basisOrEcp string, basiselems []string, basis, defstring string) string {
	if basiselems == nil { //no atoms to add basis to, do nothing
		return defstring
	}
	for _, elem := range basiselems {
		defstring = fmt.Sprintf("%s%s \"%s\" %s\n", defstring, basisOrEcp, strings.ToLower(elem), basis)
	}
	return defstring
}

//modifies the coord file such as to freeze the atoms in the slice frozen.
func (O *TMHandle) addFrozen(frozen []int) error {
	f, err := os.Open("coord")
	if err != nil {
		return Error{noCosmoPrep, Turbomole, O.inputname, err.Error(), []string{"os.Open", "addFrozen"}, true}

	}
	lines := make([]string, 0, 200) //200 is just a guess for the number of lines in the coord file
	c := bufio.NewReader(f)
	for err == nil {
		var line string
		line, err = c.ReadString('\n')
		lines = append(lines, line)
	}
	f.Close() //I cant defer it because I need it closed now.
	out, err := os.Create("coord")
	defer out.Close()
	for key, i := range lines {
		if isInInt(frozen, key-1) {
			j := strings.Replace(i, "\n", " f\n", -1)
			if _, err := fmt.Fprintf(out, j); err != nil {
				return Error{noCosmoPrep, Turbomole, O.inputname, err.Error(), []string{"fmt.Fprintf", "addFrozen"}, true}

			}
		} else {
			if _, err := fmt.Fprintf(out, i); err != nil {
				return Error{noCosmoPrep, Turbomole, O.inputname, err.Error(), []string{"fmt.Fprintf", "addFrozen"}, true}

			}
		}
	}
	return nil
}



func copy2pipe(pipe io.ReadCloser, file *os.File, end chan bool) {
	io.Copy(file, pipe)
	end <- true
}

//BuildInput builds an input for TM based int the data in atoms, coords and C.
//returns only error.
//Note that at this point the interface does not support multiplicities different from 1 and 2.
//The number in atoms is simply ignored.
func (O *TMHandle) BuildInput(coords *v3.Matrix, atoms chem.AtomMultiCharger, Q *Calc) error {
	const noDefine = "goChem/QM: Unable to run define"
	const nox2t = "goChem/QM: Unable to run x2t"
	err := os.Mkdir(O.inputname, os.FileMode(0755))
	for i := 0; err != nil; i++ {
		if strings.Contains(err.Error(), "file exists") {
			O.inputname = fmt.Sprintf("%s%d", O.inputname, i)
			err = os.Mkdir(O.inputname, os.FileMode(0755))
		} else {
			return Error{"goChem/QM: Unable to build input", Turbomole, O.inputname, err.Error(), []string{"os.Mkdir", "BuildInput"}, true}
		}
	}
	_ = os.Chdir(O.inputname)
	defer os.Chdir("..")
	//Set the coordinates in a slightly stupid way.
	chem.XYZFileWrite("file.xyz", coords, atoms)
	x2t := exec.Command("x2t", "file.xyz")
	stdout, err := x2t.StdoutPipe()
	if err != nil {
		return Error{nox2t, Turbomole, O.inputname, err.Error(), []string{"exec.StdoutPipe", "BuildInput"}, true}
	}
	coord, err := os.Create("coord")
	if err != nil {
		return Error{nox2t, Turbomole, O.inputname, err.Error(), []string{"os.Create", "BuildInput"}, true}

	}
	if err := x2t.Start(); err != nil {
		return Error{nox2t, Turbomole, O.inputname, err.Error(), []string{"exec.Start", "BuildInput"}, true}

	}
	//	var end chan bool
	//	go copy2pipe(stdout, coord, end)
	//	<-end
	io.Copy(coord, stdout)
	coord.Close()                           //not defearable
	defstring := "\n\n\na coord\nired\n*\n" //reduntant internals
	if Q.CartesianOpt {
		defstring = "\n\n\na coord\n*\nno\n"
	}
	if atoms == nil || coords == nil {
		return Error{ErrMissingCharges, Turbomole, O.inputname, "", []string{"BuildInput"}, true}
	}
	if Q.Basis == "" {
		log.Printf("no basis set assigned for TM calculation, will used the default %s, \n", O.defbasis)
		Q.Basis = O.defbasis
	}
	defstring = defstring + "b all " + Q.Basis + "\n"
	if Q.LowBasis != "" && len(Q.LBElements) > 0 {
		defstring = O.addBasis("b", Q.LBElements, Q.LowBasis, defstring)
	}
	if Q.HighBasis != "" && len(Q.HBElements) > 0 {
		defstring = O.addBasis("b", Q.HBElements, Q.HighBasis, defstring)
	}
	//Manually adding ECPs seem to be problematic, so I don't advise to do so.
	if Q.ECP != "" && len(Q.ECPElements) > 0 {
		defstring = O.addBasis("ecp", Q.ECPElements, Q.ECP, defstring)
	}
	defstring = defstring + "\n*\n"
	//The following needs to be added because some atoms (I haven't tried so many, but
	//so far only copper) causes define to ask an additional question. If one doesn't add "y\n"
	//for each of those questions, the whole input for define will be wrong.
	stupid := ""
	stupidatoms := "rting unless there is seriously unhealthy shit going on, but this is...kinda that. you condition yourself to feel bad for wanting sex, to not pursue or engage with your partner, and to repress your sexuality. I've had experiences with adults out of relationships like that, and they have no idea how to have normal sexual interactions, from flirting to pillowtalk. It's all acting, and it damages people. That said, no one should be guilted into having sex when they don't want to. Ahem.Cu" //if you want to add more stupid atoms just add then to the string: "Cu Zn"
	for i := 0; i < atoms.Len(); i++ {
		if stupidatoms == "" {
			break
		}
		if strings.Contains(stupidatoms, atoms.Atom(i).Symbol) {
			stupidatoms = strings.Replace(stupidatoms, atoms.Atom(i).Symbol, "", -1)
			stupid = stupid + "y\n"
		}
	}
	//Here we only produce singlet and doublet states (sorry). I will most certainly *not* deal with the "joys"
	//of setting other multiplicities in define.
	defstring = fmt.Sprintf("%seht\n%sy\ny\n%d\n\n", defstring, stupid, atoms.Charge()) //I add one additional "y\n"
	method, ok := tMMethods[Q.Method]
	if !ok {
		fmt.Fprintf(os.Stderr, "no method assigned for TM calculation, will used the default %s, \n", O.defmethod)
		Q.Method = O.defmethod
		Q.RI = true
	} else {
		Q.Method = method
	}
	//We only support HF and DFT
	O.command = "dscf"
	if Q.Method != "hf" {
		grid := ""
		if Q.Grid != 0 && Q.Grid <= 7 {
			grid = fmt.Sprintf("grid\n m%d\n", Q.Grid)
		}
		defstring = defstring + "dft\non\nfunc " + Q.Method + "\n" + grid + "*\n"
		if Q.RI {
			mem := 500
			if Q.Memory != 0 {
				mem = Q.Memory
			}
			defstring = fmt.Sprintf("%sri\non\nm %d\n*\n", defstring, mem)
			O.command = "ridft"
		}
	}
	defstring = O.addMARIJ(defstring, atoms, Q)
	defstring = defstring + "*\n"
	log.Println(defstring)

	//set the frozen atoms (only cartesian constraints are supported)
	if err := O.addFrozen(Q.CConstraints); err != nil {
		return errDecorate(err, "BuildInput")
	}
	if O.dryrun {
		return nil
	}
	def := exec.Command("define")
	pipe, err := def.StdinPipe()
	if err != nil {
		return Error{noDefine, Turbomole, O.inputname, err.Error(), []string{"exec.StdinPipe", "BuildInput"}, true}
	}
	defer pipe.Close()
	pipe.Write([]byte(defstring))
	if err := def.Run(); err != nil {
		return Error{noDefine, Turbomole, O.inputname, err.Error(), []string{"exec.Run", "BuildInput"}, true}
	}
	jc := jobChoose{}
	jc.opti = func() {
		O.command = "jobex"
		if Q.RI {
			O.command = O.command + " -c 200 -ri"
		} else {
			O.command = O.command + " -c 200"
		}
	}
	jc.forces = func() {
		O.command = "NumForce"
		if Q.RI {
			O.command = O.command + " -ri"
	}
	Q.Job.Do(jc)

	//Now modify control
	args := make([]string, 1, 2)
	args[0], ok = tMDisp[Q.Dispersion]
	if !ok {
		fmt.Fprintf(os.Stderr, "Dispersion correction requested not supported, will used the default: D3, \n")
		args[0] = "$disp3"
	}
	if Q.Gimic {
		O.command = "mpshift"
		args = append(args, "$gimic")
	}
	if err := O.addToControl(args, Q); err != nil {
		return errDecorate(err, "BuildInput")
	}

	//Finally the cosmo business.
	err = O.addCosmo(Q.Dielectric)
	if err != nil {
		return errDecorate(err, "BuildInput")
	}
	return nil
}


var tMMethods = map[string]string{
	"HF":     "hf",
	"hf":     "hf",
	"b3lyp":  "b3-lyp",
	"B3LYP":  "b3-lyp",
	"b3-lyp": "b3-lyp",
	"PBE":    "pbe",
	"pbe":    "pbe",
	"TPSS":   "tpss",
	"TPSSh":  "tpssh",
	"tpss":   "tpss",
	"tpssh":  "tpssh",
	"BP86":   "b-p",
	"b-p":    "b-p",
	"blyp":   "b-lyp",
	"BLYP":   "b-lyp",
	"b-lyp":  "b-lyp",
	"b97-3c": "b97-3c",
}

var tMDisp = map[string]string{
	"":       "",
	"nodisp": "",
	"D":      "$olddisp",
	"D2":     "$disp2",
	"D3":     "$disp3",
	"D3BJ":   "$disp3 -bj",
}


//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//This is a Unix-only function.
func (O *TMHandle) Run(wait bool) err error {
	os.Chdir(O.inputname)
	defer os.Chdir("..")
	filename := strings.Fields(O.command)
	//fmt.Println("nohup " + O.command + " > " + filename[0] + ".out")
	command := exec.Command("sh", "-c", "nohup "+O.command+" >"+filename[0]+".out")
	if wait == true {
		err = command.Run()
	} else {
		err = command.Start()
	}
	if err != nil {
		err = Error{ErrNotRunning, Turbomole, O.inputname, err.Error(), []string{"exec.Run/Start", "Run"}, true}

	}
	return err
}


//Energy returns the energy from the corresponding calculation, in kcal/mol.
func (O *TMHandle) Energy() (float64, error) {
	os.Chdir(O.inputname)
	defer os.Chdir("..")
	f, err := os.Open("energy")
	if err != nil {
		return 0, Error{ErrNoEnergy, Turbomole, O.inputname, err.Error(), []string{"os.Open", "Energy"}, true}
	}
	defer f.Close()
	fio := bufio.NewReader(f)
	line, err := getSecondToLastLine(fio)
	if err != nil {
		return 0, errDecorate(err, "Energy "+O.inputname)
	}
	en := strings.Fields(line)[1]
	energy, err := strconv.ParseFloat(en, 64)
	if err != nil {
		err = Error{ErrNoEnergy, Turbomole, O.inputname, err.Error(), []string{"strconv.ParseFloat", "Energy"}, true}
	}
	return energy * chem.H2Kcal, err
}

//OptimizedGeometry returns the coordinates for the optimized structure.
func (O *TMHandle) OptimizedGeometry(atoms chem.Atomer) (*v3.Matrix, error) {
	const not2x = "unable to run t2x "
	os.Chdir(O.inputname)
	defer os.Chdir("..")
	x2t := exec.Command("t2x")
	stdout, err := x2t.StdoutPipe()
	if err != nil {
		return nil, Error{ErrNoGeometry, Turbomole, O.inputname, not2x + err.Error(), []string{"exec.StdoutPipe", "OptimizedGeometry"}, true}
	}
	if err := x2t.Start(); err != nil {
		return nil, Error{ErrNoGeometry, Turbomole, O.inputname, not2x + err.Error(), []string{"exec.Start", "OptimizedGeometry"}, true}

	}
	mol, err := chem.XYZRead(stdout)
	if err != nil {
		return nil, errDecorate(err, "qm.OptimizedGeometry "+Turbomole+" "+O.inputname)
	}
	return mol.Coords[len(mol.Coords)-1], nil

}

//Gets the second to last line in a turbomole energy file given as a bufio.Reader.
//expensive on the CPU but rather easy on the memory, as the file is read line by line.
func getSecondToLastLine(f *bufio.Reader) (string, error) {
	prevline := ""
	line := ""
	var err error
	for {
		line, err = f.ReadString('\n')
		if err != nil {
			break
		}
		if !strings.Contains(line, "$end") {
			prevline = line
		} else {
			break
		}
	}
	return prevline, Error{err.Error(), Turbomole, "", "Unknown", []string{"getSecondToLastLine"}, true}
}




