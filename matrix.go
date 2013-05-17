// +build !part
/*
 * matrix.go, part of gochem.
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

import "github.com/skelterjohn/go.matrix"
import "math"
import "fmt"

//These ones are basic math, belonging more to the go.matrix package
//If there is something similar already made
//in go.matrix this functions will be deleted. Otherwise they could be
//made methods for DenseMatrix and included in go.matrix

//Some of this functions don't return error messages because they are meant to
//Be inserted in mathematical expressions and thus they need to return only one value.






//Cross3D Takes 2 3-len column or row vectors and returns a column or a row
//vector, respectively, with the Cross product of them.
//should panic
func Cross3D(a, b *matrix.DenseMatrix) (*CoordMatrix, error) {
	ar,ac := a.Dims()
	br,bc := b.Dims()
	if ac != bc || ar != br || (ac!=3 && ar!=3){
		return nil, fmt.Errorf("Malformed vectors for cross product")
	}
	if ac != 3 {
		c:=Zeros(ac,ar)
		c.T(a)
		d:=Zeros(bc,br)
		d.T(b)
		e := Cross3DRow(c, d)
		a.T(e)   //careful here, may need testing
		return a, nil
	}
	if ar != 3 {
		return Cross3DRow(a, b), nil
	}
	panic("Unreachable")
}

//Cross3DRow returns the cross product of 2 row vectors. No error checking!
func Cross3DRow(a, b *CoordMatrix) *CoordMatrix {
	vec := make([]float64, 3, 3)
	vec[0] = a.At(0, 1)*b.At(0, 2) - a.At(0, 2)*b.At(0, 1)
	vec[1] = a.At(0, 2)*b.At(0, 0) - a.At(0, 0)*b.At(0, 2)
	vec[2] = a.At(0, 0)*b.At(0, 1) - a.At(0, 1)*b.At(0, 0)
	return NewCoord(vec, 1, 3)
}

//InvSqrt return the inverse of the square root of val, or zero if
//|val|<appzero. Returns -1 if val is negative 
func InvSqrt(val float64) float64 {
	if math.Abs(val) <= appzero { //Not sure if need the 
		return 0
	} else if val < 0 { //negative
		panic("attempted to get the square root of a negative number")
	}
	return 1.0 / math.Sqrt(val)
}

//KronekerDelta is a naive implementation of the kroneker delta function.
func KronekerDelta(a, b, epsilon float64) float64 {
	if epsilon<0{
		epsilon=appzero
		}
	if math.Abs(a-b) <= epsilon {
		return 1
	}
	return 0
}





