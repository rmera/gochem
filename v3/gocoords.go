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

package v3


import (
	"fmt"
	"math"
	"strings"

	"github.com/gonum/matrix/mat64"
)

const appzero float64 = 0.000000000001 //used to correct floating point
//errors. Everything equal or less than this is considered zero. This probably sucks.

//Returns a zero-filled Matrix with vecs vectors and 3 in the other dimension.
func Zeros(vecs int) *Matrix {
	const cols int = 3
	f := make([]float64, cols*vecs, cols*vecs)
	return &Matrix{mat64.NewDense(vecs, cols,f)}
}

//METHODS

func (F *Matrix) SwapVecs(i, j int) {
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
func (F *Matrix) AddVec(A, vec *Matrix) {
	ar, ac := A.Dims()
	rr, rc := vec.Dims()
	fr, fc := F.Dims()
	if ac != rc || rr != 1 || ac != fc || ar != fr {
		panic(mat64.ErrShape)
	}
	for i := 0; i < ar; i++ {
		j := A.VecView(i)
		f := F.VecView(i)
		f.Add(j, vec)
	}
}

func (F *Matrix) DelRow(A *Matrix, i int) {
	F.DelVec(A, i)
}

func (F *Matrix) DelVec(A *Matrix, i int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if i >= ar || fc != ac || fr != (ar-1) {
		panic(mat64.ErrShape)
	}
	tempA1 := A.View(0, 0, i, ac)
	tempF1 := F.View(0, 0, i, ac)
	tempF1.Copy(tempA1)
	//now the other part
	//	if i != ar-1 {
	fmt.Println("options",ar,i,ar-i-1)
	if i<ar-1{
		tempA2 := A.View(i+1, 0, ar-i-1, ac) //The magic happens here
		tempF2 := F.View(i, 0, ar-i-1, fc)
		tempF2.Copy(tempA2)
		}
}

//return the number of vecs in F.
func (F *Matrix) NVecs() int {   //NOTE Probably just "Vecs" is a better name
	r, c := F.Dims()
	if c != 3 {
		panic(not3xXMatrix)
	}
	return r

}

//Scale each coordinates in A by the coordinate in coord.
//The result is put in F.
func (F *Matrix) ScaleByVec(A, coord *Matrix) {
	F.ScaleByVec(A, coord)
}

//Set the vectors whith index n = each value on clist, in the received to the
//n vector of A.
func (F *Matrix) SetVecs(A *Matrix, clist []int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if ac != fc || fr < len(clist) || ar < len(clist) {
		panic(mat64.ErrShape)
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
func (F *Matrix) SomeVecs(A *Matrix, clist []int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if ac != fc || fr != len(clist) || ar < len(clist) {
		panic(mat64.ErrShape)
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
func (F *Matrix) SomeVecsSafe(A *Matrix, clist []int) (error) {
	var err error
	defer func() {
		if r := recover(); r != nil {
			switch e := r.(type) {
			case Error:
				err = e
			case mat64.Error:
				err =  Error(fmt.Sprintf("%goChem/v3: Error in a gonum function: %s",e))
			default:
				panic(r)
			}
		}
	}()
	F.SomeVecs(A, clist)
	return err
}

//puts in F a matrix consistent of A over B or A to the left of B.
func (F *Matrix) StackVec(A, B *Matrix) {
	//NOTE DELETION CANDIDATE DELCAN
	F.Stack(A, B)
}

//Returns a neat string representation of a Matrix
func (F *Matrix) String() string {
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
//work if A and row reference to the same Matrix.
func (F *Matrix) SubVec(A, vec *Matrix) {
	vec.Scale(-1, vec)
	F.AddVec(A, vec)
	vec.Scale(-1, vec)
}

//Cross puts the cross product of the first vecs of a and b in the first vec of F. Panics if error.
func (F *Matrix) Cross(a, b *Matrix) {
	if a.NVecs() < 1 || b.NVecs() < 1 || F.NVecs() < 1 {
		panic("Invalid  Matrix!")
	}
	//I ask for Matrix instead of Matrix, even though  I only need the At method.
	//This is so I dont need to ensure that the rows are taken, and thus I dont need to break the
	//API if the matrices become col-major.
	F.Set(0, 0, a.At(0, 1)*b.At(0, 2)-a.At(0, 2)*b.At(0, 1))
	F.Set(0, 1, a.At(0, 2)*b.At(0, 0)-a.At(0, 0)*b.At(0, 2))
	F.Set(0, 2, a.At(0, 0)*b.At(0, 1)-a.At(0, 1)*b.At(0, 0))
}

//METHODS Not Vec specific.

//AddFloat puts in the receiver a matrix which elements are
//those of matrix A plus the float B.
func (F *Matrix) AddFloat(A *Matrix, B float64) {
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
//reference to the same Matrix.
func (F *Matrix) AddRow(A, row *Matrix) {
	F.AddVec(A, row)
}

//ScaleByCol scales each column of matrix A by Col, putting the result
//in the received.
func (F *Matrix) ScaleByCol(A, Col mat64.Matrix) {
	ar, ac := A.Dims()
	cr, cc := Col.Dims()
	fr, fc := F.Dims()
	if ar != cr || cc > 1 || ar != fr || ac != fc {
		panic(mat64.ErrShape)
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
func (F *Matrix) ScaleByRow(A, Row *Matrix) {     //NOTE it should be called ScaleByVec
	ar, ac := A.Dims()
	rr, rc := Row.Dims()
	fr, fc := F.Dims()
	if ac != rc || rr != 1 || ar != fr || ac != fc {
		panic(mat64.ErrShape)
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
func (F *Matrix) RowView(i int) *Matrix {
	return F.VecView(i)
}

//SubRow subtracts the row vector row to each row of the matrix A, putting
//the result on the receiver. Panics if matrices are mismatched.  It will not
//work if A and row reference to the same Matrix.
func (F *Matrix) SubRow(A, row *Matrix) {
	F.SubVec(A, row)
}

func (F *Matrix) Unit(A *Matrix) {
	if A.Dense != F.Dense {
		F.Copy(A)
	}
	norm := 1.0 / F.Norm(0)
	F.Scale(norm, F)
}



//KronekerDelta is a naive implementation of the kroneker delta function.
func KronekerDelta(a, b, epsilon float64) float64 {
	if epsilon < 0 {
		epsilon = appzero
	}
	if math.Abs(a-b) <= epsilon {
		return 1
	}
	return 0
}
