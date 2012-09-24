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

	
package chem

import 


type IntConstraint struct{
	Kind byte
	Atoms []int
	}

type PointCharge{
	Charge float64
	Coords *matrix.DenseMatrix
	}

type QMCalc struct {
	Method string
	Basis string
	AuxBasis string //for RI calculations
	AuxColBasis string //for RICOSX or similar calculations
	HighBasis string  //a bigger basis for certain atoms
	Lowbasis string  //lower basis for certain atoms
	HBatoms []int
	LBatoms []int
	CConstraint []int //cartesian contraints
	IConstraints []IntCostraint //internal constraints
	Dielectric float64
	Solventmethod string
	Disperssion string  //D, D2, D3
	Others string //analysis methods, etc
	PCharges []PointCharge
	Optimize bool
	NCPU int
	Execute bool //Execute (true), or just produce the inputs (false)?
	Command string //this is the whole command to execute
	Guess string //initial guess
	OldMO bool //Try to look for a file with MO. The 
	
	}


func BuildOrcaInput(atoms Reference, coords matrix.DenseMatrix, C QMCalc, path string) error{
	if c.Method=="" || basis==""{
		return fmt.Errorf("Not enough options for the optimization")
		}
	disp=""
	switch C.Disperssion{
		case "":
		case "D2":
		disp="VDW"
		case "D3":
		disp="VDW3"
		default:
		disp="VDW3"
		}
	if atoms==nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
		}
	opt:=""
	if C.Optimise==true{
		opt=="Opt"
		}
	if C.OldMO==true{
		//set something to look for a *.gbw in the path
		C.Guess="MORead"
		//set some variable to 
		}
	MainOptions:=[]string{"!",C.Method,C.basis,C.AuxBasis,C.AuxColBasis,C.Guess,opt,C.Others}
	mainline:=strings.Join(MainOptions," ")
	}
	



















import "fmt"

func main() {

}
