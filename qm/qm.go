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
 * Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche
 *
 */

package qm

import (
	"fmt"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

// builds an input for a QM calculation
type InputBuilder interface {
	//Sets the name for the job, used for input
	//and output files. The extentions will depend on the program.
	SetName(name string)

	//BuildInput builds an input for the QM program based int the data in
	//atoms, coords and C. returns only error.
	BuildInput(coords *v3.Matrix, atoms chem.AtomMultiCharger, Q *Calc) error
}

// Runs a QM calculation
type Runner interface {
	//Run runs the QM program for a calculation previously set.
	//it waits or not for the result depending of the value of
	//wait.
	Run(wait bool) (err error)
}

// Builds inputs and runs a QM calculations
type BuilderRunner interface {
	InputBuilder
	Runner
}

// Allows to recover energy and optimized geometries from a QM calculation
type EnergyGeo interface {

	//Energy gets the last energy for a  calculation by parsing the
	//QM program's output file. Return error if fail. Also returns
	//Error ("Probable problem in calculation")
	//if there is a energy but the calculation didnt end properly.
	Energy() (float64, error)

	//OptimizedGeometry reads the optimized geometry from a calculation
	//output. Returns error if fail. Returns Error ("Probable problem
	//in calculation") if there is a geometry but the calculation didnt
	//end properly*
	OptimizedGeometry(atoms chem.Atomer) (*v3.Matrix, error)
}

// Handle is an interface for a mostly-full functionality QM program
// where "functionality" reflects it's degree of support in goChem
type Handle interface {
	BuilderRunner
	EnergyGeo
}

const (
	ErrProbableProblem = "goChem/QM: Probable problem with calculations" //this is never to be used for fatal errors
	ErrMissingCharges  = "goChem/QM: Missing charges or coordinates"
	ErrCantValue       = "goChem/QM: Can't obtain requested value" //for whatever doesn't match anything else
	ErrNoEnergy        = "goChem/QM: No energy in output"
	ErrNoFreeEnergy    = "goChem/QM: No free energy in output. Forces calculation might have not been performed"
	ErrNoCharges       = "goChem/QM: Unable to read charges from  output"
	ErrNoGeometry      = "gochem/QM: Unable to read geometry from output"
	ErrNotRunning      = "gochem/QM: Couldn't run calculation"
	ErrCantInput       = "goChem/QM: Can't build input file"
)

const (
	Orca      = "Orca"
	Mopac     = "Mopac"
	Turbomole = "Turbomole"
	NWChem    = "NWChem"
	Fermions  = "Fermions++"
	XTB       = "XTB" //this may go away if Orca starts supporting XTB.
)

//errors

// Error represents a decorable QM error.
type Error struct {
	message    string
	code       string //the name of the QM program giving the problem, or empty string if none
	inputname  string //the input file that has problems, or empty string if none.
	additional string
	deco       []string
	critical   bool
}

// Error returns a string with an error message.
func (err Error) Error() string {
	return fmt.Sprintf("%s (%s/%s) Message: %s", err.message, err.inputname, err.code, err.additional)
}

// Decorate will add the dec string to the decoration slice of strings of the error,
// and return the resulting slice.
func (err Error) Decorate(dec string) []string {
	err.deco = append(err.deco, dec)
	return err.deco
}

// Code returns the name of the program that ran/was meant to run the
// calculation that caused the error.
func (err Error) Code() string { return err.code } //May not be needed

// InputName returns the name of the input file which processing caused the error
func (err Error) InputName() string { return err.inputname }

// Critical return whether the error is critical or it can be ifnored
func (err Error) Critical() bool { return err.critical }

// errDecorate is a helper function that asserts that the error is
// implements chem.Error and decorates the error with the caller's name before returning it.
// if used with a non-chem.Error error, it will cause a panic.
func errDecorate(err error, caller string) error {
	err2 := err.(chem.Error) //I know that is the type returned byt initRead
	err2.Decorate(caller)
	return err2
}

//end errors

// jobChoose is a structure where each QM handler has to provide a closure that makes the proper arrangements for each supported case.
type jobChoose struct {
	opti    func()
	forces  func()
	sp      func()
	md      func()
	charges func()
}

// Job is a structure that define a type of calculations.
// The user should set one of these to true,
// and goChem will see that the proper actions are taken. If the user sets more than one of the
// fields to true, the priority will be Opti>Forces>SP (i.e. if Forces and SP are true,
// only the function handling forces will be called).
type Job struct {
	Opti    bool
	Forces  bool
	SP      bool
	MD      bool
	Charges bool
}

// Do sets the job set to true in J, according to the corresponding function in plan. A "nil" plan
// means that the corresponding job is not supported by the QM handle and we will default to single point.
func (J *Job) Do(plan jobChoose) {
	if J == nil {
		return
	}
	//now the actual options
	if J.Opti {
		plan.opti()
		return
	}
	if J.Forces && plan.forces != nil {
		plan.forces()
		return
	}
	if J.MD && plan.md != nil {
		plan.md()
		return
	}
	if J.Charges && plan.charges != nil { //the default option is a single-point
		plan.charges()
		return
	}

	if plan.sp != nil { //the default option is a single-point
		plan.sp()
		return
	}

}

//Container for constraints to internal coordinates
//type IntConstraint struct {
//	Kind  byte
//	Atoms []int
//}

// PointCharge is a container for a point charge, such as those used in QM/MM
// calculations
type PointCharge struct {
	Charge float64
	Coords *v3.Matrix
}

// IConstraint is a container for a constraint to internal coordinates
type IConstraint struct {
	CAtoms []int
	Val    float64
	Class  byte // B: distance, A: angle, D: Dihedral
	UseVal bool //if false, don't add any value to the constraint (which should leave it at the value in the starting structure. This migth not work on every program, but it works in ORCA.
}

// Calc is a structure for the general representation of a calculation
// mostly independent of the QM program (although, of course, some methods will not work in some programs)
type Calc struct {
	Method       string
	Basis        string
	RI           bool
	RIJ          bool
	CartesianOpt bool   //Do the optimization in cartesian coordinates.
	BSSE         string //Correction for BSSE
	auxBasis     string //for RI calculations
	auxColBasis  string //for RICOSX or similar calculations
	HighBasis    string //a bigger basis for certain atoms
	LowBasis     string //lower basis for certain atoms
	HBAtoms      []int
	LBAtoms      []int
	HBElements   []string
	LBElements   []string
	CConstraints []int //cartesian contraints
	IConstraints []*IConstraint
	ECPElements  []string //list of elements with ECP.
	//	IConstraints []IntConstraint //internal constraints
	Dielectric float64
	//	Solventmethod string
	Dispersion string //D2, D3, etc.
	Others     string //analysis methods, etc
	//	PCharges []PointCharge
	Guess string //initial guess
	Grid  int
	OldMO bool //Try to look for a file with MO.
	Job   *Job //NOTE: This should probably be a pointer: FIX!  NOTE2: Fixed it, but must check and fix whatever is now broken.
	//The following 3 are only for MD simulations, will be ignored in every other case.
	MDTime       int     //simulation time (whatever unit the program uses!)
	MDTemp       float64 //simulation temperature (K)
	MDPressure   int     //simulation pressure (whatever unit the program uses!)
	SCFTightness int     //1,2 and 3. 0 means the default
	OptTightness int
	SCFConvHelp  int
	ECP          string //The ECP to be used. It is the programmers responsibility to use a supported ECP (for instance, trying to use 10-electron core ECP for Carbon will fail)
	Gimic        bool
	NBO          bool
	Memory       int //Max memory to be used in MB (the effect depends on the QM program)
}

// Utilities here
func (Q *Calc) SetDefaults() {
	Q.RI = true
	//	Q.BSSE = "gcp"
	Q.Dispersion = "D3"
	Q.Memory = 3000
}

// isIn is a helper for the RamaList function,
// returns true if test is in container, false otherwise.
func isInInt(container []int, test int) bool {
	if container == nil {
		return false
	}
	for _, i := range container {
		if test == i {
			return true
		}
	}
	return false
}

// Same as the previous, but with strings.
func isInString(container []string, test string) bool {
	if container == nil {
		return false
	}
	for _, i := range container {
		if test == i {
			return true
		}
	}
	return false
}
