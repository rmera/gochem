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
	}


func OrcaOptimizer(Reference, )




















import "fmt"

func main() {

}
