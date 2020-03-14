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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 *
 */
//In order to use this part of the library you need the xtb program, which must be obtained from Prof. Stefan Grimme's group.
//Please cite the the xtb references if you used the program.

/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package qm

import (
	//	"bufio"
	"fmt"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
)

//Note that the default methods and basis vary with each program, and even
//for a given program they are NOT considered part of the API, so they can always change.
type XTBHandle struct {
	command   string
	inputname string
	nCPU      int
	options   []string
}

func NewXTBHandle() *XTBHandle {
	run := new(XTBHandle)
	run.SetDefaults()
	return run
}

//XTBHandle methods

//Sets the number of CPU to be used
func (O *XTBHandle) SetnCPU(cpu int) {
	O.nCPU = cpu
}

func (O *XTBHandle) Command() string {
	return O.command
}

func (O *XTBHandle) SetName(name string) {
	O.inputname = name
}

func (O *XTBHandle) SetCommand(name string) {
	O.command = name
}

func (O *XTBHandle) SetDefaults() {
	O.command = os.ExpandEnv("xtb")
	//	if O.command == "/xtb" { //if XTBHOME was not defined
	//		O.command = "./xtb"
	//	}
	cpu := runtime.NumCPU() / 2
	O.nCPU = cpu

}

//BuildInput builds an input for XTB. Right now it's very limited, only singlets are allowed and
//only unconstrained optimizations and single-points.
func (O *XTBHandle) BuildInput(coords *v3.Matrix, atoms chem.AtomMultiCharger, Q *Calc) error {

	//Only error so far
	if atoms == nil || coords == nil {
		return Error{ErrMissingCharges, "XTB", O.inputname, "", []string{"BuildInput"}, true}
	}
	err := chem.XYZFileWrite(O.inputname+".xyz", coords, atoms)
	if err != nil {
		return Error{ErrCantInput, "XTB", O.inputname, "", []string{"BuildInput"}, true}
	}
	//	mem := ""
	if Q.Memory != 0 {
		//Here we can adjust memory if needed
	}
	//Now lets write the thing
	if O.inputname == "" {
		O.inputname = "gochem"
	}
	O.options = make([]string, 2, 6)
	O.options[0] = O.command
	O.options[1] = O.inputname + ".xyz"
	O.options = append(O.options, fmt.Sprintf("-c %d", atoms.Charge()))
	O.options = append(O.options, fmt.Sprintf("-u %d", (atoms.Multi()-1)))
	if O.nCPU > 1 {
		O.options = append(O.options, fmt.Sprintf("-P %d", O.nCPU))
	}

	if Q.Dielectric > 0 {
		solvent, ok := dielectric2Solvent[int(Q.Dielectric)]
		if ok {
			O.options = append(O.options, "-g "+solvent)
		}
	}
	//O.options = append(O.options, "-gfn")
	fixed := ""
	if Q.CConstraints != nil {
		fixed = "atoms: "
		for _, v := range Q.CConstraints {
			fixed = fixed + strconv.Itoa(v) + ", " //0-based indexes
		}
		strings.TrimRight(fixed, ",")
		fixed = fixed + "\n"
		xcontrol, err := os.Create("xcontrol")
		if err != nil {
			return err
		}
		xcontrol.Write([]byte("$fix\n"))
		xcontrol.Write([]byte("force constant=10000\n"))
		xcontrol.Write([]byte(fixed))
		xcontrol.Close()
		O.options = append(O.options, "-I xcontrol")
	}
	jc := jobChoose{}
	jc.opti = func() {
		O.options = append(O.options, "-o normal")
	}

	Q.Job.Do(jc)
	return nil

}

//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//Not waiting for results works
//only for unix-compatible systems, as it uses bash and nohup.
func (O *XTBHandle) Run(wait bool) (err error) {
	if wait == true {
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
		command := exec.Command("sh", "-c", O.command+fmt.Sprintf(" %s.xyz %s > %s.out  2>&1", O.inputname, strings.Join(O.options[2:], " "), O.inputname))
		err = command.Run()

	} else {
		command := exec.Command("sh", "-c", "nohup "+O.command+fmt.Sprintf(" %s.xyz %s > %s.out &", O.inputname, O.options[1:], O.inputname))
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

//Reads the latest geometry from an XTB optimization. It doesn't actually need the chem.Atomer
//but requires it so XTBHandle fits with the QM interface.
func (O *XTBHandle) OptimizedGeometry(atoms chem.Atomer) (*v3.Matrix, error) {
	if !O.normalTermination() {
		return nil, Error{ErrNoGeometry, XTB, O.inputname, "Calculation didn't end normally", []string{"OptimizedGeometry"}, true}
	}
	mol, err := chem.XYZFileRead("xtbopt.xyz") //Trying to run several calculations in parallel in the same directory will fail as the output has always the same name.
	if err != nil {
		return nil, Error{ErrNoGeometry, XTB, O.inputname, "", []string{"OptimizedGeometry"}, true}
	}
	return mol.Coords[0], nil
}

//This checks that an xtb calculation has terminated normally
//I know this duplicates code, I wrote this one first and then the other one.
func (O *XTBHandle) normalTermination() bool {
	if searchBackwards("normal termination of x", fmt.Sprintf("%s.out", O.inputname)) != "" || searchBackwards("abnormal termination of x", fmt.Sprintf("%s.out", O.inputname)) == "" {
		return true
	}
	//	fmt.Println(fmt.Sprintf("%s.out", O.inputname), searchBackwards("normal termination of x",fmt.Sprintf("%s.out", O.inputname))) ////////////////////
	return false
}

//search a file backwards, i.e., starting from the end, for a string. Returns the line that contains the string, or an empty string.
func searchBackwards(str, filename string) string {
	var ini int64 = 0
	var end int64 = 0
	var first bool
	first = true
	buf := make([]byte, 1)
	//	fmt.Println("no wei", filename) ////////////////////////
	f, err := os.Open(filename)
	if err != nil {
		//		fmt.Println(err.Error())	 ////////////////
		return ""
	}
	defer f.Close()
	var i int64 = 1
	for ; ; i++ {
		if _, err := f.Seek(-1*i, 2); err != nil {
			//		fmt.Println(err.Error()) ///////////////////
			return ""
		}
		if _, err := f.Read(buf); err != nil {
			//		fmt.Println(err.Error()) //////////////////
			return ""
		}
		if buf[0] == byte('\n') && first == false {
			first = true
		} else if buf[0] == byte('\n') && end == 0 {
			end = i
		} else if buf[0] == byte('\n') && ini == 0 {
			ini = i
			f.Seek(-1*(ini), 2)
			bufF := make([]byte, ini-end)
			f.Read(bufF)
			//		fmt.Println("vieja", string(bufF))////////////////////
			if strings.Contains(string(bufF), str) {
				return string(bufF)
			}
			//	first=false
			end = 0
			ini = 0
		}

	}
	return ""
}

//Gets the energy of a previous XTB calculations.
//Returns error if problem, and also if the energy returned that is product of an
//abnormally-terminated ORCA calculation. (in this case error is "Probable problem
//in calculation")
func (O *XTBHandle) Energy() (float64, error) {
	var err error
	var energy float64
	energyline := searchBackwards("total E       :", fmt.Sprintf("%s.out", O.inputname))
	if energyline == "" {
		return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"searchBackwards", "Energy"}, true}
	}
	split := strings.Fields(energyline)
	if len(split) < 4 {
		return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"Energy"}, true}

	}
	energy, err = strconv.ParseFloat(split[3], 64)
	if err != nil {
		return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"strconv.ParseFloat", "Energy"}, true}
	}

	return energy * chem.H2Kcal, err //dummy thin
}

var dielectric2Solvent = map[int]string{
	80: "h2o",
	5:  "chcl3",
	9:  "ch2cl2",
	21: "acetone",
	37: "acetonitrile",
	33: "methanol",
	2:  "toluene",
	7:  "thf",
	47: "dmso",
	38: "dmf",
}
