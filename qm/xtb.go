/*
 * xtb.go, part of gochem.
 *
 *
 * Copyright 2016 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 *
 */
//In order to use this part of the library you need the xtb program, which must be obtained from Prof. Stefan Grimme's group.
//Please cite the the xtb references if you used the program.

package qm

import (
	//	"bufio"
	"bufio"
	"fmt"
	"math"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

// XTBHandle represents an xtb calculation
type XTBHandle struct {
	//Note that the default methods and basis vary with each program, and even
	//for a given program they are NOT considered part of the API, so they can always change.
	//This is unavoidable, as methods change with time
	command        string
	inputname      string
	nCPU           int
	options        []string
	gfnff          bool
	relconstraints bool
	force          float64
	wrkdir         string
	inputfile      string
}

// NewXTBHandle initializes and returns an xtb handle
// with values set to their defaults. Defaults might change
// as new methods appear, so they are not part of the API.
func NewXTBHandle() *XTBHandle {
	run := new(XTBHandle)
	run.SetDefaults()
	return run
}

//XTBHandle methods

// SetnCPU sets the number of CPU to be used
func (O *XTBHandle) SetnCPU(cpu int) {
	O.nCPU = cpu
}

// Command returns the path and name for the xtb excecutable
func (O *XTBHandle) Command() string {
	return O.command
}

// SetName sets the name for the calculations
// which is defines the input and output file names
func (O *XTBHandle) SetName(name string) {
	O.inputname = name
}

// SetCommand sets the path and name for the xtb excecutable
func (O *XTBHandle) SetCommand(name string) {
	O.command = name
}

// SetWorkDir sets the name of the working directory for the calculations
func (O *XTBHandle) SetWorkDir(d string) {
	O.wrkdir = d
}

// RelConstraints sets the use of relative contraints
// instead of absolute position restraints
// with the force force constant. If force is
// less than 0, the default value is employed.
func (O *XTBHandle) RelConstraints(force float64) {
	if force > 0 {
		//goChem units (Kcal/A^2) to xtb units (Eh/Bohr^2)
		O.force = force * (chem.Kcal2H / (chem.A2Bohr * chem.A2Bohr))
	}
	O.relconstraints = true

}

// SetDefaults sets calculations parameters to their defaults.
// Defaults might change
// as new methods appear, so they are not part of the API.
func (O *XTBHandle) SetDefaults() {
	O.command = os.ExpandEnv("xtb")
	//	if O.command == "/xtb" { //if XTBHOME was not defined
	//		O.command = "./xtb"
	//	}
	cpu := runtime.NumCPU() / 2
	O.nCPU = cpu

}

func (O *XTBHandle) seticonstraints(Q *Calc, xcontrol []string) []string {
	//Here xtb expects distances in A and angles in deg, so no conversion needed.
	var g2x = map[byte]string{
		'B': "distance: ",
		'A': "angle: ",
		'D': "dihedral: ",
	}

	for _, v := range Q.IConstraints {
		constra := g2x[v.Class]

		for _, w := range v.CAtoms {
			constra += strconv.Itoa(w+1) + ", " //1-based indexes for xtb
		}
		strings.TrimRight(constra, ",")
		if v.UseVal {
			constra += fmt.Sprintf(" %4.2f\n", v.Val)
		} else {
			constra += " auto\n"
		}
		xcontrol = append(xcontrol, constra)

	}
	return xcontrol
}

// BuildInput builds an input for XTB. Right now it's very limited, only singlets are allowed and
// only unconstrained optimizations and single-points.
func (O *XTBHandle) BuildInput(coords *v3.Matrix, atoms chem.AtomMultiCharger, Q *Calc) error {
	//Now lets write the thing
	if O.wrkdir != "" && !strings.HasSuffix(O.wrkdir, "/") {
		O.wrkdir += "/"
	}
	w := O.wrkdir
	if O.inputname == "" {
		O.inputname = "gochem"
	}
	//Only error so far
	if atoms == nil || coords == nil {
		return Error{ErrMissingCharges, "XTB", O.inputname, "", []string{"BuildInput"}, true}
	}
	err := chem.XYZFileWrite(w+O.inputname+".xyz", coords, atoms)
	if err != nil {
		return Error{ErrCantInput, "XTB", O.inputname, "", []string{"BuildInput"}, true}
	}
	//	mem := ""
	if Q.Memory != 0 {
		//Here we can adjust memory if needed
	}

	xcontroltxt := make([]string, 0, 10)
	O.options = make([]string, 0, 6)
	O.options = append(O.options, O.command)
	if Q.Method == "gfnff" {
		O.gfnff = true
	}
	O.options = append(O.options, O.inputname+".xyz")
	O.options = append(O.options, fmt.Sprintf("-c %d", atoms.Charge()))
	O.options = append(O.options, fmt.Sprintf("-u %d", (atoms.Multi()-1)))
	if O.nCPU > 1 {
		O.options = append(O.options, fmt.Sprintf("-P %d", O.nCPU))
	}
	//Added new things to select a method in xtb
	if !isInString([]string{"gfn1", "gfn2", "gfn0", "gfnff"}, Q.Method) {
		O.options = append(O.options, "--gfn 2") //default method
	} else if Q.Method != "gfnff" {
		m := strings.ReplaceAll(Q.Method, "gfn", "") //so m should be "0", "1" or "2"
		O.options = append(O.options, "--gfn "+m)    //default method
	}

	if Q.Dielectric > 0 && Q.Method != "gfn0" { //as of the current version, gfn0 doesn't support implicit solvation
		solvent, ok := dielectric2Solvent[int(Q.Dielectric)]
		if ok {
			O.options = append(O.options, "--alpb "+solvent)
		}
	}
	//O.options = append(O.options, "-gfn")
	fixed := ""
	fixtoken := "$fix\n"
	if O.relconstraints {
		force := ""
		if O.force > 0 {
			force = fmt.Sprintf("force constant = %4.2f\n", O.force)
		}
		fixtoken = "$constrain\n" + force
	}
	if Q.CConstraints != nil || Q.IConstraints != nil {
		xcontroltxt = append(xcontroltxt, fixtoken)
		if Q.CConstraints != nil {
			fixed = "atoms: "
			for _, v := range Q.CConstraints {
				fixed = fixed + strconv.Itoa(v+1) + ", " //1-based indexes
			}
			strings.TrimRight(fixed, ",")
			fixed = fixed + "\n"
			xcontroltxt = append(xcontroltxt, fixed)
		}
		if Q.IConstraints != nil {
			if !O.relconstraints {
				xcontroltxt = append(xcontroltxt, ("$end\n$constrain\n"))
			}
			xcontroltxt = O.seticonstraints(Q, xcontroltxt)

		}
		xcontroltxt = append(xcontroltxt, "$end\n")

	}
	jc := jobChoose{}
	jc.opti = func() {
		add := "-o normal"
		if Q.OptTightness == 2 {
			add = "-o tight"
		} else if Q.OptTightness > 2 {
			add = "-o verytight"
		}
		O.options = append(O.options, add)
	}
	jc.forces = func() {
		O.options = append(O.options, "--ohess")
	}

	jc.md = func() {
		O.options = append(O.options, "--md")
		//There are specific settings needed with gfnff, mainly, a shorter timestep
		//The restart=false option doesn't have any effect, but it's added so it's easier later to use sed or whatever to change it to true, and  restart
		//a calculation.
		if Q.Method == "gfnff" {
			xcontroltxt = append(xcontroltxt, fmt.Sprintf("$md\n temp=%5.3f\n time=%d\n velo=false\n nvt=true\n step=2.0\n hmass=4.0\n shake=0\n restart=false\n$end", Q.MDTemp, Q.MDTime))
		} else {
			xcontroltxt = append(xcontroltxt, fmt.Sprintf("$md\n temp=%5.3f\n time=%d\n velo=false\n nvt=true\n restart=false\n$end", Q.MDTemp, Q.MDTime))
		}

	}
	//	O.options = append(O.options, "--input xcontrol")
	O.options = append(O.options, Q.Others)
	Q.Job.Do(jc)
	if len(xcontroltxt) == 0 {
		return nil //no need to write a control file
	}
	O.inputfile = O.inputname + ".inp" //if not input file was written
	//this will just be an empty string.
	xcontrol, err := os.Create(w + O.inputfile)
	if err != nil {
		return err
	}
	for _, v := range xcontroltxt {
		xcontrol.WriteString(v)

	}
	xcontrol.Close()
	return nil
}

// Run runs the command given by the string O.command
// it waits or not for the result depending on wait.
// Not waiting for results works
// only for unix-compatible systems, as it uses bash and nohup.
func (O *XTBHandle) Run(wait bool) (err error) {
	var com string
	extraoptions := ""
	if len(O.options) >= 3 {
		extraoptions = strings.Join(O.options[2:], " ")
	}
	inputfile := ""
	if O.inputfile != "" {
		inputfile = fmt.Sprintf("--input %s", O.inputfile)
	}

	if O.gfnff {
		com = fmt.Sprintf(" --gfnff %s.xyz   %s  %s > %s.out  2>&1", O.inputname, inputfile, extraoptions, O.inputname)
	} else {

		com = fmt.Sprintf(" %s.xyz  %s  %s > %s.out  2>&1", O.inputname, inputfile, extraoptions, O.inputname)
	}
	if wait {
		//It would be nice to have this logging as an option.
		//log.Printf(O.command + com) //this is stderr, I suppose
		command := exec.Command("sh", "-c", O.command+com)
		command.Dir = O.wrkdir
		err = command.Run()

	} else {
		command := exec.Command("sh", "-c", "nohup "+O.command+com)
		command.Dir = O.wrkdir
		err = command.Start()
	}
	if err != nil {
		err = Error{ErrNotRunning, XTB, O.inputname, err.Error(), []string{"exec.Start", "Run"}, true}
	}
	if err != nil {
		return err
	}
	os.Remove("xtbrestart")
	return nil
}

// OptimizedGeometry returns the latest geometry from an XTB optimization. It doesn't actually need the chem.Atomer
// but requires it so XTBHandle fits with the QM interface.
func (O *XTBHandle) OptimizedGeometry(atoms chem.Atomer) (*v3.Matrix, error) {
	inp := O.wrkdir + O.inputname
	if !O.normalTermination() {
		return nil, Error{ErrNoGeometry, XTB, inp, "Calculation didn't end normally", []string{"OptimizedGeometry"}, true}
	}
	mol, err := chem.XYZFileRead(O.wrkdir + "xtbopt.xyz") //Trying to run several calculations in parallel in the same directory will fail as the output has always the same name.
	if err != nil {
		return nil, Error{ErrNoGeometry, XTB, inp, "", []string{"OptimizedGeometry"}, true}
	}
	return mol.Coords[0], nil
}

// This checks that an xtb calculation has terminated normally
// I know this duplicates code, I wrote this one first and then the other one.
func (O *XTBHandle) normalTermination() bool {
	inp := O.wrkdir + O.inputname
	if searchBackwards("normal termination of x", fmt.Sprintf("%s.out", inp)) != "" || searchBackwards("abnormal termination of x", fmt.Sprintf("%s.out", inp)) == "" {
		return true
	}
	return false
}

// search a file backwards, i.e., starting from the end, for a string. Returns the line that contains the string, or an empty string.
// I really really should have commented this one.
func searchBackwards(str, filename string) string {
	var ini int64 = 0
	var end int64 = 0
	var first bool
	first = true
	buf := make([]byte, 1)
	f, err := os.Open(filename)
	if err != nil {
		return ""
	}
	defer f.Close()
	var i int64 = 1
	for ; ; i++ {
		if _, err := f.Seek(-1*i, 2); err != nil {
			return ""
		}
		if _, err := f.Read(buf); err != nil {
			return ""
		}
		if buf[0] == byte('\n') && first == false {
			first = true
		} else if buf[0] == byte('\n') && end == 0 {
			end = i
		} else if buf[0] == byte('\n') && ini == 0 {
			i--
			ini = i
			f.Seek(-1*(ini), 2)
			bufF := make([]byte, ini-end)
			f.Read(bufF)
			if strings.Contains(string(bufF), str) {
				return string(bufF)
			}
			//	first=false
			end = 0
			ini = 0
		}

	}
}

// Energy returns the energy of a previous XTB calculations, in kcal/mol.
// Returns error if problem, and also if the energy returned that is product of an
// abnormally-terminated ORCA calculation. (in this case error is "Probable problem
// in calculation")
func (O *XTBHandle) Energy() (float64, error) {
	inp := O.wrkdir + O.inputname
	var err error
	var energy float64
	energyline := searchBackwards("TOTAL ENERGY", fmt.Sprintf("%s.out", inp))
	if energyline == "" {
		return 0, Error{ErrNoEnergy, XTB, inp, fmt.Sprintf("%s.out", inp), []string{"searchBackwards", "Energy"}, true}
	}
	split := strings.Fields(energyline)
	if len(split) < 5 {
		return 0, Error{ErrNoEnergy, XTB, inp, err.Error(), []string{"Energy"}, true}

	}
	energy, err = strconv.ParseFloat(split[3], 64)
	if err != nil {
		return 0, Error{ErrNoEnergy, XTB, inp, err.Error(), []string{"strconv.ParseFloat", "Energy"}, true}
	}

	return energy * chem.H2Kcal, err //dummy thin
}

// LargestImaginary returns the absolute value of the wave number (in 1/cm) for the largest imaginary mode in the vibspectrum file
// produced by a forces calculation with xtb. Returns an error and -1 if unable to check.
func (O *XTBHandle) LargestImaginary() (float64, error) {
	largestimag := 0.0
	vibf, err := os.Open(O.wrkdir + "vibspectrum")
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return -1, Error{ErrCantValue, XTB, "vibspectrum", err.Error(), []string{"os.Open", "LargestImaginary"}, true}
	}
	vib := bufio.NewReader(vibf)
	for i := 0; i < 3; i++ {
		_, err := vib.ReadString('\n') //The text in "data" could be anything, including just "\n"
		if err != nil {
			return -1, Error{ErrCantValue, XTB, "vibspectrum", err.Error(), []string{"ReadString", "LargestImaginary"}, true}

		}
	}
	for {
		line, err := vib.ReadString('\n')
		if err != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			if strings.Contains(err.Error(), "EOF") {
				err = nil //it's not an actual error
				break
			} else {
				return -1, Error{ErrCantValue, XTB, "vibspectrum", err.Error(), []string{"ReadString", "LargestImaginary"}, true}

			}
		}

		fields := strings.Fields(line)
		if len(fields) < 5 {
			return -1, Error{ErrCantValue, XTB, "vibspectrum", "Can't parse vibspectrum", []string{"ReadString", "LargestImaginary"}, true}
		}
		wave, err := strconv.ParseFloat(fields[len(fields)-4], 64)
		if err != nil {
			return -1, Error{ErrCantValue, XTB, "vibspectrum", "Can't parse vibspectrum", []string{"strconv.ParseFloat", "LargestImaginary"}, true}
		}
		if wave > 0.0 {
			return largestimag, nil //no more imaginary frequencies so we just return the largest so far.
		} else if math.Abs(wave) > largestimag {
			largestimag = math.Abs(wave)
		}
	}
	return largestimag, nil
}

// FixImaginary prepares and runs a calculation on a geometry, produced by xtb on a previous Hessian calculation, which
// is distorted along the main imaginary mode found, if any. It such mode was not found, and thus the geometry was not
// produced by xtb, FixImaginary returns an error.
func (O *XTBHandle) FixImaginary(wait bool) error {
	var com string
	var err error
	if _, err := os.Stat("xtbhess.coord"); os.IsNotExist(err) {
		return fmt.Errorf("xtbhess.coord doesn't exist. There is likely no significant imaginary mode")
	}
	if O.gfnff {
		com = fmt.Sprintf(" --gfnff xtbhess.coord  --input %s.inp  %s > %s.out  2>&1", O.inputname, strings.Join(O.options[2:], " "), O.inputname)
	} else {

		com = fmt.Sprintf(" xtbhess.coord  --input %s.inp  %s > %s.out  2>&1", O.inputname, strings.Join(O.options[2:], " "), O.inputname)
	}
	if wait == true {
		//log.Printf(com) //this is stderr, I suppose
		command := exec.Command("sh", "-c", O.command+com)
		command.Dir = O.wrkdir
		err = command.Run()

	} else {
		command := exec.Command("sh", "-c", "nohup "+O.command+com)
		command.Dir = O.wrkdir
		err = command.Start()
	}
	if err != nil {
		err = Error{ErrNotRunning, XTB, O.inputname, err.Error(), []string{"exec.Start", "Run"}, true}
	}
	if err != nil {
		return err
	}
	os.Remove("xtbrestart")
	return nil

}

// FreeEnergy returns the Gibbs free energy of a previous XTB calculations.
// A frequencies/solvation calculation is needed for this to work. FreeEnergy does _not_ check that the structure was at a minimum. You can check that with
// the LargestIm
func (O *XTBHandle) FreeEnergy() (float64, error) {
	var err error
	var energy float64
	inp := O.wrkdir + O.inputname
	energyline := searchBackwards("total free energy", fmt.Sprintf("%s.out", inp))
	if energyline == "" {
		return 0, Error{ErrNoFreeEnergy, XTB, inp, fmt.Sprintf("%s.out", inp), []string{"searchBackwards", "FreeEnergy"}, true}
	}
	split := strings.Fields(energyline)
	if len(split) < 4 {
		return 0, Error{ErrNoFreeEnergy, XTB, inp, err.Error(), []string{"Energy"}, true}

	}
	energy, err = strconv.ParseFloat(split[4], 64)
	if err != nil {
		return 0, Error{ErrNoFreeEnergy, XTB, inp, err.Error(), []string{"strconv.ParseFloat", "Energy"}, true}
	}

	return energy * chem.H2Kcal, err //err should be nil at this point.
}

// MDAverageEnergy gets the average potential and kinetic energy along a trajectory.
func (O *XTBHandle) MDAverageEnergy(start, skip int) (float64, float64, error) {
	inp := O.wrkdir + O.inputname
	var potential, kinetic float64
	if !O.normalTermination() {
		return 0, 0, Error{ErrNoEnergy, XTB, inp, "Calculation didn't end normally", []string{"MDAverageEnergy"}, true}
	}
	outname := fmt.Sprintf("%s.out", inp)
	outfile, err := os.Open(outname)
	if err != nil {
		return 0, 0, Error{ErrNoEnergy, XTB, inp, "Couldn't open output file", []string{"MDAverageEnergy"}, true}
	}
	out := bufio.NewReader(outfile)
	reading := false
	cont := 0
	read := 0
	for {
		line, err := out.ReadString('\n')
		//	fmt.Println("LINE", line) /////////
		if err != nil && strings.Contains(err.Error(), "EOF") {
			break
		} else if err != nil {
			return 0, 0, Error{ErrNoEnergy, XTB, inp, "Error while iterating through output file", []string{"MDAverageEnergy"}, true}
		}
		if strings.Contains(line, "time (ps)    <Epot>      Ekin   <T>   T     Etot") {
			reading = true
			continue
		}
		if !reading {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) != 7 {
			continue
		}
		cont++
		if (cont-1)%skip != 0 || (cont-1) < start {
			continue
		}
		K, err := strconv.ParseFloat(fields[3], 64)
		if err != nil {
			return 0, 0, Error{ErrNoEnergy, XTB, inp, fmt.Sprintf("Error while retrieving %d th kinetic energy", cont), []string{"MDAverageEnergy"}, true}

		}
		V, err := strconv.ParseFloat(fields[3], 64)
		if err != nil {
			return 0, 0, Error{ErrNoEnergy, XTB, inp, fmt.Sprintf("Error while retrieving %d th potential energy", cont), []string{"MDAverageEnergy"}, true}

		}
		fmt.Println("potential", V) //////////
		kinetic += K
		potential += V
		read++
	}
	N := float64(read)
	if math.IsNaN(potential/N) || math.IsNaN(kinetic/N) { //note that we still return whatever we got here, in addition to the error. The user can decide.
		return potential / N, kinetic / N, Error{ErrProbableProblem, XTB, inp, "At least one of the energies is NaN", []string{"MDAverageEnergy"}, true}
	}
	return potential / N, kinetic / N, nil
}

var dielectric2Solvent = map[int]string{
	80:   "h2o",
	5:    "chcl3",
	9:    "ch2cl2",
	10:   "octanol",
	21:   "acetone",
	37:   "acetonitrile",
	33:   "methanol",
	2:    "toluene",
	1:    "hexadecane", //not quite 1
	7:    "thf",
	47:   "dmso",
	38:   "dmf",
	1000: "woctanol", //This is a hackish way to have both dry and wet octanol. I gave wet octanol an fake epsilon that won't be used by anything else.
	//really, what I should do is to add to the API a way to specify either epsilon or solvent name. FIX
}

//old code

//		out, err := os.Create(fmt.Sprintf("%s.out", O.inputname))
//		if err != nil {
//			return Error{ErrNotRunning, XTB, O.inputname, "", []string{"Run"}, true}
//		}
//		ferr, err := os.Create(fmt.Sprintf("%s.err", O.inputname))
//
//		if err != nil {
//			return Error{ErrNotRunning, XTB, O.inputname, "", []string{"Run"}, true}
//		}
//		defer out.Close()
//		defer ferr.Close()
//		fullCommand:=strings.Join(O.options," ")
//		fmt.Println(fullCommand) //("Command", O.command, O.options) ////////////////////////
//		command := exec.Command(fullCommand) //, O.options...)
//		command.Stdout = out
//		command.Stderr = ferr
//		err = command.Run()
//		fmt.Println(O.command+fmt.Sprintf(" %s.xyz %s > %s.out &", O.inputname, strings.Join(O.options[2:]," "), O.inputname)) ////////////////////////
