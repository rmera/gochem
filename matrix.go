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
func Cross3D(a, b *matrix.DenseMatrix) (*matrix.DenseMatrix, error) {
	ac := a.Cols()
	ar := a.Rows()
	bc := b.Cols()
	br := b.Rows()
	if ac != bc || ar != br {
		return nil, fmt.Errorf("Malformed vectors for cross product")
	}
	if ac != 3 {
		//Ok, Im sure one can do this better.
		c := a.Transpose()
		d := b.Transpose()
		e := Cross3DRow(c, d)
		f := e.Transpose()
		return f, nil
	}
	if ar != 3 {
		return Cross3DRow(a, b), nil
	}
	panic("Unreachable")
}

//Cross3DRow returns the cross product of 2 row vectors. No error checking!
func Cross3DRow(a, b *matrix.DenseMatrix) *matrix.DenseMatrix {
	vec := make([]float64, 3, 3)
	vec[0] = a.Get(0, 1)*b.Get(0, 2) - a.Get(0, 2)*b.Get(0, 1)
	vec[1] = a.Get(0, 2)*b.Get(0, 0) - a.Get(0, 0)*b.Get(0, 2)
	vec[2] = a.Get(0, 0)*b.Get(0, 1) - a.Get(0, 1)*b.Get(0, 0)
	return matrix.MakeDenseMatrix(vec, 1, 3)
}

//InvSqrt return the inverse of the square root of val, or zero if
//|val|<appzero. Returns -1 if val is negative 
func InvSqrt(val float64) float64 {
	if math.Abs(val) <= appzero { //Not sure if need the 
		return 0
	} else if val < 0 { //negative
		return -1 //might change
	}
	return 1.0 / math.Sqrt(val)
}

//KronekerDelta is a naive implementation of the kroneker delta function.
func KronekerDelta(a, b float64) float64 {
	if math.Abs(a-b) <= appzero {
		return 1
	}
	return 0
}




//DMScaleByCol scales each column of matrix A by Col.
func DMScaleByCol(A, Col *matrix.DenseMatrix) error {
	Rows := A.Rows()
	if Rows != Col.Rows() || Col.Cols() > 1 {
		return fmt.Errorf("DMScaleByCol: Malformed matrices for scaling")
	}
	for i := 0; i < Rows; i++ {
		A.ScaleRow(i, Col.Get(i, 0))
	}
	return nil
}

//DMScaleByRow each row of matrix A by Row.
func DMScaleByRow(A, Row *matrix.DenseMatrix) error {
	Cols := A.Cols()
	if Cols != Row.Cols() || Row.Rows() > 1 {
		return fmt.Errorf("DMScaleByRow: Malformed matrices for scaling")
	}
	for i := 0; i < Cols; i++ {
		mult := Row.Get(0, i)
		for j := 0; j < A.Rows(); j++ {
			A.Set(j, i, mult*A.Get(j, i))
		}
	}
	return nil
}


//DMaddfloat returns a matrix which elements are those of matrix A plus the float B.
func DMaddfloat(A *matrix.DenseMatrix, B float64) *matrix.DenseMatrix {
	Rows := A.Rows()
	Cols := A.Cols()
	copy := matrix.MakeDenseCopy(A)
	for i := 0; i < Cols; i++ {
		for j := 0; j < Rows; j++ {
			copy.Set(j, i, (A.Get(j, i) + B))
		}
	}
	return copy
}

//DMDelRow removes Row i from matrix A. 
func DMDelRow(A *matrix.DenseMatrix, i int) (*matrix.DenseMatrix, error) {
	//Im not sure if copying takes place.
	var err error
	lenght := A.Rows()
	upper := A.GetMatrix(0, 0, i, 3)
	lower := A.GetMatrix(i+1, 0, lenght-i-1, 3)
	A, err = upper.Stack(lower)
	return A, err
}
