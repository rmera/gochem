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
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
	"bufio"
	"github.com/rmera/gochem"
	"github.com/rmera/gochem/v3"
)

//Note that the default methods and basis vary with each program, and even
//for a given program they are NOT considered part of the API, so they can always change.
type XTBHandle struct {
	command     string
	inputname   string
	nCPU        int
	options		string
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

func (O *XTBHandle) SetName(name string) {
	O.inputname = name
}

func (O *XTBHandle) SetCommand(name string) {
	O.command = name
}


func (O *XTBHandle) SetDefaults() {
	O.command = os.ExpandEnv("${XTBHOME}/xtb")
	if O.command == "/xtb" { //if ORCA_PATH was not defined
		O.command = "./xtb"
	}
	cpu := runtime.NumCPU()
	O.nCPU = cpu

}

//BuildInput builds an input for XTB. Right now it's very limited, only singlets are allowed and 
//only unconstrained optimizations and single-points.
func (O *XTBHandle) BuildInput(coords *v3.Matrix, atoms chem.AtomMultiCharger, Q *Calc) error {
	//Only error so far
	if atoms == nil || coords == nil {
		return Error{ErrMissingCharges, "XTB", O.inputname, "", []string{"BuildInput"}, true}
		}
		err:=chem.XYZFileWrite(O.inputname+".xyz",coords,atoms)
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
	f, err := os.OpenFile(O.inputname+".xyz", os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
    	return Error{ErrCantInput, "XTB", O.inputname, "", []string{"BuildInput"}, true}
	}

	if _, err = f.WriteString(fmt.Sprintf("$set\ncub_cal 1\nchrg %2d\n$end",atoms.Charge())); err != nil {
    	return Error{ErrCantInput, "XTB", O.inputname, "", []string{"BuildInput"}, true} //it would be nice to differenciate this error from the previous.
	}

	f.Close() //Won't use defer, as we need this file written and saved before this function exits.
	jc := jobChoose{}
	jc.opti = func() {
		O.options="-opt"
	}
	jc.sp = func() {
		O.options="-sp"
	}

	Q.Job.Do(jc)
	if Q.Dielectric > 0 {
		O.options=O.options+" -gbsa h2o" //Only water supported for now
	}

	return nil
}

//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//Not waiting for results works
//only for unix-compatible systems, as it uses bash and nohup.
func (O *XTBHandle) Run(wait bool) (err error) {
	if wait == true {
		out, err := os.Create(fmt.Sprintf("%s.out", O.inputname))
		if err != nil {
			return Error{ErrNotRunning, XTB, O.inputname, "", []string{"Run"}, true}
		}
		defer out.Close()
		command := exec.Command(O.command, fmt.Sprintf("%s.xyz %s",O.inputname,O.options))
		command.Stdout = out
		err = command.Run()

	} else {
		command := exec.Command("sh", "-c", "nohup "+O.command+fmt.Sprintf(" %s.xyz %s > %s.out &", O.inputname, O.options, O.inputname))
		err = command.Start()
	}
	if err != nil {
		err = Error{ErrNotRunning, XTB, O.inputname, err.Error(), []string{"exec.Start", "Run"}, true}
	}
	return err
}

//Reads the latest geometry from an XTB optimization. It doesn't actually need the chem.Atomer
//but requires it so XTBHandle fits with the QM interface.
func (O *XTBHandle) OptimizedGeometry(atoms chem.Atomer) (*v3.Matrix, error) {
	mol,err:=chem.XYZFileRead("xtbopt.coord") //Trying to run several calculations in parallel in the same directory will fail as the output has always the same name.
	return mol.Coords[0], err
}

//Gets the energy of a previous XTB calculations.
//Returns error if problem, and also if the energy returned that is product of an
//abnormally-terminated ORCA calculation. (in this case error is "Probable problem
//in calculation")
func (O *XTBHandle) Energy() (float64, error) {
	var err error
	var energy float64
	file, err := os.Open(fmt.Sprintf("energy"))
	if err != nil {
		return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"os.Open", "Energy"}, true}
	}
	defer file.Close()
	out := bufio.NewReader(file)
	err = Error{ErrNoEnergy, XTB, O.inputname, "", []string{"Energy"}, true}
	var line string
	line, err = out.ReadString('\n')
		if err != nil {
			return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"os.file.ReadString", "Energy"}, true}
		}
	line, err = out.ReadString('\n')
		if err != nil {
			return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"os.file.ReadString", "Energy"}, true}
		}
	split := strings.Fields(line)
	if len(split) < 4 {
		err = Error{ErrNoEnergy, Mopac, O.inputname, "", []string{"Energy"}, true}
			return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"Energy"}, true}

		}
	energy, err = strconv.ParseFloat(split[1], 64)
	if err != nil {
		return 0, Error{ErrNoEnergy, XTB, O.inputname, err.Error(), []string{"strconv.ParseFloat", "Energy"}, true}
	}



	return energy * chem.H2Kcal, err //dummy thin
}

