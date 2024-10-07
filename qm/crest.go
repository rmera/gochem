/*
 * crest.go, part of gochem.
 *
 *
 * Copyright 2024 Raul Mera <rmeraa{at}academicosdotutadotcl
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
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

// CrestHandle represents an xtb calculation
type CrestHandle struct {
	//Note that the default methods and basis vary with each program, and even
	//for a given program they are NOT considered part of the API, so they can always change.
	//This is unavoidable, as methods change with time
	command        string
	inputname      string
	nCPU           int
	options        []string
	relconstraints bool
	wrkdir         string
	inputfile      string
	RunType        string     //entropy, protonate, deprotonate, search (default)
	Temperatures   [3]float64 //initial, final, step
	EThres         float64
	RMSDThres      float64
}

// NewCrestHandle initializes and returns an xtb handle
// with values set to their defaults. Defaults might change
// as new methods appear, so they are not part of the API.
func NewCrestHandle() *CrestHandle {
	run := new(CrestHandle)
	run.SetDefaults()
	return run
}

//CrestHandle methods

// SetnCPU sets the number of CPU to be used
func (O *CrestHandle) SetnCPU(cpu int) {
	O.nCPU = cpu
}

// Command returns the path and name for the xtb excecutable
func (O *CrestHandle) Command() string {
	return O.command
}

// SetName sets the name for the calculations
// which is defines the input and output file names
func (O *CrestHandle) SetName(name string) {
	O.inputname = name
}

// SetCommand sets the path and name for the xtb excecutable
func (O *CrestHandle) SetCommand(name string) {
	O.command = name
}

// SetWorkDir sets the name of the working directory for the calculations
func (O *CrestHandle) SetWorkDir(d string) {
	O.wrkdir = d
}

// SetDefaults sets calculations parameters to their defaults.
// Defaults might change
// as new methods appear, so they are not part of the API.
func (O *CrestHandle) SetDefaults() {
	O.command = os.ExpandEnv("crest")
	//	if O.command == "/xtb" { //if CrestHOME was not defined
	//		O.command = "./xtb"
	//	}
	cpu := runtime.NumCPU() / 2
	O.nCPU = cpu

}

// BuildInput builds an input for Crest. Right now it's very limited, only singlets are allowed and
// only unconstrained optimizations and single-points.
func (O *CrestHandle) BuildInput(coords *v3.Matrix, atoms chem.AtomMultiCharger, Q *Calc) error {
	errid := "CrestHandle/BuildInput"
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
		return fmt.Errorf("%s: no molecule or coordinates given", errid)
	}
	err := chem.XYZFileWrite(w+O.inputname+".xyz", coords, atoms)
	if err != nil {
		return fmt.Errorf("%s: Couldn't write xyz file: %w ", errid, err)
	}
	//	mem := ""

	O.options = make([]string, 0, 6)
	//	O.options = append(O.options, O.command)
	O.options = append(O.options, fmt.Sprintf("--chrg %d", atoms.Charge()))
	O.options = append(O.options, fmt.Sprintf("--uhf %d", (atoms.Multi()-1)))
	if O.nCPU > 1 {
		O.options = append(O.options, fmt.Sprintf("-P %d", O.nCPU))
	}
	//Added new things to select a method in xtb
	if !isInString([]string{"gfn1", "gfn2", "gfn0", "gfnff"}, Q.Method) {
		O.options = append(O.options, "--gfn2") //default method
	} else {
		O.options = append(O.options, "--"+Q.Method) //default method
	}
	if Q.Dielectric > 0 && Q.Method != "gfn0" { //as of the current version, gfn0 doesn't support implicit solvation
		solvent, ok := dielectric2Solvent[int(Q.Dielectric)]
		if ok {
			O.options = append(O.options, "--alpb "+solvent)
		}
	}

	o := "--optlev vtight"
	if Q.OptTightness > 0 {
		if Q.OptTightness < 2 {
			o = "--optlev normal"
		}
		if Q.OptTightness == 2 {
			o = "--optlev tight"
		}
	}
	O.options = append(O.options, o)

	//run type. It's your responsibility to set only one of them to true.

	switch strings.ToLower(O.RunType) {
	case "entropy":
		O.options = append(O.options, "--entropy")
	case "protonate":
		O.options = append(O.options, "--protonate")
	case "deprotonate":
		O.options = append(O.options, "--deprotonate")
	case "v2":
		O.options = append(O.options, "--v2")
	case "v1":
		O.options = append(O.options, "--v1")
	case "v3":
		O.options = append(O.options, "--v3")
	case "v4":
		O.options = append(O.options, "--v4")
	case "tautomerize":
		O.options = append(O.options, "--tautomerize")
	case "":
		O.options = append(O.options, "") //I could have just left this empty
	default:
		log.Printf("%s: Runtype %s not available, will do the CREST default (currently, v3)", errid, O.RunType)
	}

	ts := O.Temperatures
	if ts == [3]float64{0, 0, 0} {
		ts = [3]float64{298.15, 299.15, 1}
	}
	O.options = append(O.options, fmt.Sprintf("--trange %5.2f %5.2f %5.2f", ts[0], ts[1], ts[2]))

	if O.EThres > 0 { //crest expect this options in kcal, so no conversion needed
		O.options = append(O.options, fmt.Sprintf("--ewin %4.1f", O.EThres))
	}
	if O.RMSDThres > 0 { //crest expect this options in kcal, so no conversion needed
		O.options = append(O.options, fmt.Sprintf("--rthr %4.1f", O.RMSDThres))
	}

	O.options = append(O.options, fmt.Sprintf("--temp %5.2f", ts[0]))

	if Q.CConstraints != nil || Q.IConstraints != nil {
		xtbh := NewXTBHandle()
		constraintsfile := "constraints"
		xtbh.SetName(constraintsfile)
		xtbh.SetWorkDir(O.wrkdir)
		err = xtbh.BuildInput(coords, atoms, Q)
		if err != nil {
			return fmt.Errorf("%s: Couldn't produce the constraints file: %w", errid, err)
		}
		O.options = append(O.options, "--cinp "+constraintsfile+".inp")
	}
	O.options = append(O.options, O.inputname+".xyz")

	return nil
}

// Run runs the command given by the string O.command
// it waits or not for the result depending on wait.
// Not waiting for results works
// only for unix-compatible systems, as it uses bash and nohup.
func (O *CrestHandle) Run(wait bool) (err error) {
	errid := "CrestHandle/Run"
	var com string
	options := ""
	if len(O.options) >= 3 {
		options = strings.Join(O.options, " ")
	}

	com = fmt.Sprintf(" %s > %s.out  2>&1", options, O.inputname)

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
		err = fmt.Errorf("%s: %w", errid, err)
	}
	os.Remove("xtbrestart")
	return nil
}

// Checks that an CREST calculation has terminated normally
func (O *CrestHandle) normalTermination() bool {
	inp := O.wrkdir + O.inputname
	if searchBackwards("CREST terminated normally", fmt.Sprintf("%s.out", inp)) != "" {
		return true
	}
	return false
}

func (O *CrestHandle) ConformerEnergies() ([]float64, error) {
	ei := "CrestHanle/ConformerEnergies"
	if !O.normalTermination() {
		return nil, fmt.Errorf("%s: CREST run didn't finish normally", ei)
	}

	finp, err := os.Open(O.wrkdir + "crest_conformers.xyz")
	if err != nil {
		return nil, fmt.Errorf("%s: %w", ei, err)
	}
	defer finp.Close()
	energies := make([]float64, 0, 10)
	fr := bufio.NewReader(finp)
	var line string
	trim := strings.TrimSpace
	var reade bool
	nenergy := 0
	for line, err = fr.ReadString('\n'); err == nil; line, err = fr.ReadString('\n') {
		if _, err2 := strconv.Atoi(trim(line)); err2 == nil {
			//I'm only trying to detect the lines after which comes the energy. Those lines
			//should have only an integer, and they should be the only lines to be like that
			reade = true
			continue
		}
		if reade {
			e, err := strconv.ParseFloat(trim(line), 64)
			if err != nil {
				return nil, fmt.Errorf("%s: Couldn't parse energy %d: %w", ei, nenergy, err)
			}
			energies = append(energies, e*chem.H2Kcal)
			nenergy++
			reade = false
			continue
		}
	}
	if errors.Is(err, io.EOF) {
		err = nil
	}
	return energies, err
}

// Returns conformers from CREST as a trajectory, returning also the first one as a molecule.
// if asmolecule is given and true, then it returns the whole trajectory in the molecule, and nil
// for the trajectory.
func (O *CrestHandle) Conformers(asmolecule ...bool) (*chem.Molecule, *chem.XYZTraj, error) {
	ei := "CrestHandle/Conformers"
	if !O.normalTermination() {
		return nil, nil, fmt.Errorf("%s: CREST run didn't finish normally", ei)
	}
	if len(asmolecule) > 0 && asmolecule[0] {
		mol, err := chem.XYZFileRead(O.wrkdir + "crest_conformers.xyz")
		if err != nil {
			err = fmt.Errorf("%s: Failed to retrieve conformers: %w", ei, err)
		}
		return mol, nil, err
	}
	mol, traj, err := chem.XYZFileAsTraj(O.wrkdir + "crest_conformers.xyz")
	if err != nil {
		err = fmt.Errorf("%s: Failed to retrieve conformers: %w", ei, err)
	}
	return mol, traj, err
}

// Returns the conformational and vibrational entropies at the first temperature given in the
// options (298.15 by default).
func (O *CrestHandle) Entropy() (float64, float64, error) {
	ei := "CrestHanle/Entropy"
	if !O.normalTermination() {
		return 0, 0, fmt.Errorf("%s: CREST run didn't finish normally", ei)
	}
	inp := O.wrkdir + O.inputname
	var err error
	var sconf float64
	var svib float64
	energyline := searchBackwards("Sconf   =", fmt.Sprintf("%s.out", inp))
	if energyline == "" {
		return 0, 0, fmt.Errorf("%s: Couldn't find conformational entropy in output", ei)
	}
	split := strings.Fields(energyline)
	if len(split) < 3 {
		return 0, 0, fmt.Errorf("%s: Bad format in found conformational entropy line: %s", ei, energyline)

	}
	sconf, err = strconv.ParseFloat(split[2], 64)
	if err != nil {
		return 0, 0, fmt.Errorf("%s: Couldn't parse conformational entropy to a number: %s %s %w", ei, energyline, split[2], err)

	}
	sline := searchBackwards("+ δSrrho  =", fmt.Sprintf("%s.out", inp))
	if sline == "" {
		return 0, 0, fmt.Errorf("%s: Couldn't find vibrational entropy in output", ei)
	}
	split = strings.Fields(sline)
	if len(split) < 4 {
		return 0, 0, fmt.Errorf("%s:  Bad format in found  vibrational entropy: %s", ei, sline)

	}
	svib, err = strconv.ParseFloat(split[3], 64)
	if err != nil {
		return 0, 0, fmt.Errorf("%s: Couldn't parse vibrational entropy to a number: %s %s %w", ei, energyline, split[3], err)
	}

	return sconf / 1000.0, svib / 1000.0, err //It appears that the entropies are given in cal/mol*K, we give then
	//in gochem units, kcal/mol*K
}
