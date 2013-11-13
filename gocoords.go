/*
 * gocoords.go, part of gochem.
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

//Package chem provides atom and molecule structures, facilities for reading and writing some
//files used in computational chemistry and some functions for geometric manipulations and shape
//indicators.

package chem

import (
	"fmt"
	"math"
	"strings"
)

// Matrix is the basic matrix interface type. This is redundant and exactly equivalent to the implementation in gonum
type Matrix interface {
	// Dims returns the dimensions of a Matrix.
	Dims() (r, c int)

	// At returns the value of a matrix element at (r, c). It will panic if r or c are
	// out of bounds for the matrix.
	At(r, c int) float64
}

const appzero float64 = 0.000000000001 //used to correct floating point
//errors. Everything equal or less than this is considered zero.

//Returns a zero-filled VecMatrix with cos vectors and 3 in the other dimension.
func ZeroVecs(cos int) *VecMatrix {
	const cols int = 3
	dens := gnZeros(cos, cols)
	return &VecMatrix{dens.Dense}
}

func chemDense2VecMatrix(A *chemDense) *VecMatrix {
	return &VecMatrix{A.Dense}
}

func vecMatrix2chemDense(A *VecMatrix) *chemDense {
	return &chemDense{A.Dense}
}

//METHODS

func (F *VecMatrix) SwapVecs(i, j int) {
	if i >= F.NVecs() || j >= F.NVecs() {
		panic("Indexes out of range")
	}
	rowi := F.Row(nil, i)
	rowj := F.Row(nil, j)
	for k := 0; k < 3; k++ {
		F.Set(i, k, rowj[k])
		F.Set(j, k, rowi[k])
	}
}

//Adds a vector to the  coordmatrix A putting the result on the received.
//depending on whether the underlying matrix to coordmatrix
//is col or row major, it could add a col or a row vector.
func (F *VecMatrix) AddVec(A, vec *VecMatrix) {
	ar, ac := A.Dims()
	rr, rc := vec.Dims()
	fr, fc := F.Dims()
	if ac != rc || rr != 1 || ac != fc || ar != fr {
		panic(gnErrShape)
	}
	for i := 0; i < ar; i++ {
		j := A.VecView(i)
		f := F.VecView(i)
		f.Add(j, vec)
	}
}

//puts a copy of matrix A without the vector i
//in the received. Vector could be a col or row vector depending
//on the majorship of the implementation.
func (F *VecMatrix) DelVec(A *VecMatrix, i int) {
	F.DelRow(A, i)
}

func (F *VecMatrix) DelRow(A *VecMatrix, i int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if i >= ar || fc != ac || fr != (ar-1) {
		panic(gnErrShape)
	}
	tempA1 := A.View(0, 0, i, ac)
	tempF1 := F.View(0, 0, i, ac)
	tempF1.Copy(tempA1)
	//now the other part
	//	if i != ar-1 {
	tempA2 := A.View(i+1, 0, ar-i-1, ac) //The magic happens here
	tempF2 := F.View(i, 0, ar-i-1, fc)
	tempF2.Copy(tempA2)
	//	}
}

//return the number of vecs in F. Panics if the
//other dimmension is not 3.
func (F *VecMatrix) NVecs() int {
	r, c := F.Dims()
	if c != 3 {
		panic(not3xXMatrix)
	}
	return r

}

//Scale each coordinates in A by the coordinate in coord.
//The result is put in F.
func (F *VecMatrix) ScaleByVec(A, coord *VecMatrix) {
	F.ScaleByVec(A, coord)
}

//Set the vectors whith index n = each value on clist, in the received to the
//n vector of A.
func (F *VecMatrix) SetVecs(A *VecMatrix, clist []int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if ac != fc || fr < len(clist) || ar < len(clist) {
		panic(gnErrShape)
	}
	for key, val := range clist {
		for j := 0; j < ac; j++ {
			F.Set(val, j, A.At(key, j))
		}
	}
}

//Returns a matrix contaning all the ith rows of matrix A,
//where i are the numbers in clist. The rows are in the same order
//than the clist.
func (F *VecMatrix) SomeVecs(A *VecMatrix, clist []int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if ac != fc || fr != len(clist) || ar < len(clist) {
		panic(gnErrShape)
	}
	for key, val := range clist {
		for j := 0; j < ac; j++ {
			F.Set(key, j, A.At(val, j))
		}
	}
}

//Returns a matrix contaning all the ith vectors of matrix A,
//where i are the numbers in clist. The vectors are in the same order
//than the clist. Returns an error instead of panicking.
func (F *VecMatrix) SomeVecsSafe(A *VecMatrix, clist []int) (err error) {
	f := func() { F.SomeVecs(A, clist) }
	return gnMaybe(gnPanicker(f))
}

//puts in F a matrix consistent of A over B or A to the left of B.
//DELCAN
func (F *VecMatrix) StackVec(A, B *VecMatrix) {
	F.Stack(A, B)
}

//Returns a neat string representation of a VecMatrix
func (F *VecMatrix) String() string {
	r, c := F.Dims()
	v := make([]string, r+2, r+2)
	v[0] = "\n["
	v[len(v)-1] = " ]"
	row := make([]float64, c, c)
	for i := 0; i < r; i++ {
		F.Row(row, i) //now row has a slice witht he row i
		if i == 0 {
			v[i+1] = fmt.Sprintf("%6.2f %6.2f %6.2f\n", row[0], row[1], row[2])
			continue
		} else if i == r-1 {
			v[i+1] = fmt.Sprintf(" %6.2f %6.2f %6.2f", row[0], row[1], row[2])
			continue
		} else {
			v[i+1] = fmt.Sprintf(" %6.2f %6.2f %6.2f\n", row[0], row[1], row[2])
		}
	}
	v[len(v)-2] = strings.Replace(v[len(v)-2], "\n", "", 1)
	return strings.Join(v, "")
}

//SubVec subtracts the vector  to each vector of the matrix A, putting
//the result on the receiver. Panics if matrices are mismatched.  It will not
//work if A and row reference to the same VecMatrix.
func (F *VecMatrix) SubVec(A, vec *VecMatrix) {
	vec.Scale(-1, vec)
	F.AddVec(A, vec)
	vec.Scale(-1, vec)
}

//Cross puts the cross product of the first vecs of a and b in the first vec of F. Panics if error.
func (F *VecMatrix) Cross(a, b *VecMatrix) {
	if a.NVecs() < 1 || b.NVecs() < 1 || F.NVecs() < 1 {
		panic("Invalid  VecMatrix!")
	}
	//I ask for VecMatrix instead of Matrix, even though  I only need the At method.
	//This is so I dont need to ensure that the rows are taken, and thus I dont need to break the
	//API if the matrices become col-major.
	F.Set(0, 0, a.At(0, 1)*b.At(0, 2)-a.At(0, 2)*b.At(0, 1))
	F.Set(0, 1, a.At(0, 2)*b.At(0, 0)-a.At(0, 0)*b.At(0, 2))
	F.Set(0, 2, a.At(0, 0)*b.At(0, 1)-a.At(0, 1)*b.At(0, 0))
}

//METHODS Not Vec specific.

//AddFloat puts in the receiver a matrix which elements are
//those of matrix A plus the float B.
func (F *VecMatrix) AddFloat(A *VecMatrix, B float64) {
	ar, ac := A.Dims()
	if F != A {
		F.Copy(A)
	}
	for i := 0; i < ar; i++ {
		for j := 0; j < ac; j++ {
			F.Set(i, j, A.At(i, j)+B)
		}
	}
}

//AddRow adds the row vector row to each row of the matrix A, putting
//the result on the receiver. Panics if matrices are mismatched. It will not work if A and row
//reference to the same VecMatrix.
func (F *VecMatrix) AddRow(A, row *VecMatrix) {
	F.AddVec(A, row)
}

//Puts A**exp on the receiver. This function could probably
//be written in a concurrent way
func (F *chemDense) Pow(A Matrix, exp float64) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if ar != fr || ac != fc {
		panic(gnErrShape)
	}
	for i := 0; i < ar; i++ {
		for j := 0; j < ac; j++ {
			F.Set(i, j, math.Pow(A.At(i, j), exp))
		}

	}
}

//ScaleByCol scales each column of matrix A by Col, putting the result
//in the received.
func (F *VecMatrix) ScaleByCol(A, Col Matrix) {
	ar, ac := A.Dims()
	cr, cc := Col.Dims()
	fr, fc := F.Dims()
	if ar != cr || cc > 1 || ar != fr || ac != fc {
		panic(gnErrShape)
	}
	if F != A {
		F.Copy(A)
	}
	for i := 0; i < ac; i++ {
		temp := F.ColView(i)
		temp.MulElem(temp, Col)
	}

}

//ScaleByRow scales each column of matrix A by Col, putting the result
//in the received.
func (F *VecMatrix) ScaleByRow(A, Row *VecMatrix) {
	ar, ac := A.Dims()
	rr, rc := Row.Dims()
	fr, fc := F.Dims()
	if ac != rc || rr != 1 || ar != fr || ac != fc {
		panic(gnErrShape)
	}
	if F != A {
		F.Copy(A)
	}
	for i := 0; i < ac; i++ {
		temp := F.RowView(i)
		temp.MulElem(temp, Row)
	}
}

//Puts a view of the given row of the matrix in the receiver
func (F *VecMatrix) RowView(i int) *VecMatrix {
	return F.VecView(i)
}

//SubRow subtracts the row vector row to each row of the matrix A, putting
//the result on the receiver. Panics if matrices are mismatched.  It will not
//work if A and row reference to the same VecMatrix.
func (F *VecMatrix) SubRow(A, row *VecMatrix) {
	F.SubVec(A, row)
}

func (F *VecMatrix) Unit(A *VecMatrix) {
	if A.Dense != F.Dense {
		F.Copy(A)
	}
	norm := 1.0 / F.Norm(0)
	F.Scale(norm, F)
}
