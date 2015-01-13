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

import (
	"github.com/rmera/gochem"
	"fmt"
	)

//builds an input for a QM calculation 
type  InputBuilder interface {
	//Sets the name for the job, used for input
	//and output files. The extentions will depend on the program.
	SetName(name string)

	//BuildInput builds an input for the QM program based int the data in
	//atoms, coords and C. returns only error.
	BuildInput(coords *chem.VecMatrix, atoms chem.ReadRef, Q *Calc) error
}


//Runs a QM calculation
type Runner interface {
	//Run runs the QM program for a calculation previously set.
	//it waits or not for the result depending of the value of
	//wait.
	Run(wait bool) (err error)

}


type BuilderRunner interface {
		InputBuilder
		Runner
}

//Allows to recover energy and optimized geometries from a QM calculation
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
	OptimizedGeometry(atoms chem.Ref) (*chem.VecMatrix, error) /*NOTE: The "Probable problem..." error should have a type, and should report the program used. Also, chem.Ref is probably not needed here. Atomer is probably enough*/

}

//This allows to set QM calculations using different programs.
type Handle interface {
	BuilderRunner
	EnergyGeo

}


const (
	ProbableProblem  = "goChem/QM: Probable problem with calculations" //this is never to be used for fatal errors
	MissingCharges   = "goChem/QM: Missing charges or coordinates"
	NoEnergy  =   "goChem/QM: No energy in output"
	NoGeometry =  "gochem/QM: Unable to read Geometry from input"
	NotRunning =  "gochem/QM: Couldn't run calculation"
	NoInput =  "goChem/QM: Can't build input file"

)

const (
	Orca = "Orca"
	Mopac = "Mopac"
	Turbomole = "Turbomole"
	NWChem = "NWChem"
	Fermions = "Fermions++"

)

type Error  struct {
	message string
	code string //the name of the QM program giving the problem, or empty string if none
	inputname string //the input file that has problems, or empty string if none.
	additional string
	critical bool
}
func (err Error) Error() string { return fmt.Sprintf("%s (%s/%s) Message: %s",err.message,err.inputname, err.code,err.additional)  }

func (err Error) Code() string {return err.code} //May not be needed

func (err Error) InputName() string {return err.inputname}







type IntConstraint struct {
	Kind  byte
	Atoms []int
}

type PointCharge struct {
	Charge float64
	Coords *chem.VecMatrix
}

type IConstraint struct {
	CAtoms []int
	Val    float64
	Class  byte // B: distance, A: angle, D: Dihedral
}

type Calc struct {
	Method       string
	Basis        string
	RI           bool
	RIJ          bool
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
	Others      string //analysis methods, etc
	//	PCharges []PointCharge
	Guess        string //initial guess
	Grid         int
	OldMO        bool //Try to look for a file with MO. The
	Optimize     bool
	SCFTightness int
	SCFConvHelp  int
	ECP          string //The ECP to be used. It is the programmers responsibility to use a supported ECP (for instance, trying to use 10-electron core ECP for Carbon will fail)
	Gimic        bool
	Memory       int //Max memory to be used in MB (the effect depends on the QM program)
}

func (Q *Calc) SetDefaults() {
	Q.RI = true
//	Q.BSSE = "gcp"
	Q.Dispersion = "D3"
}

//Utilities here

//isIn is a helper for the RamaList function,
//returns true if test is in container, false otherwise.
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

//Same as the previous, but with strings.
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
