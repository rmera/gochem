/*
 * gonum.go, part of gochem.
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
	"github.com/skelterjohn/go.matrix"
	"math"
	"sort"
)

/*Here I make a -very incomplete- implementation of the gonum api backed by go.matrix, which will enable me to port gochem to gonum.
 * Since the agreement in the gonum community was NOT to build a temporary implementation, I just build the functions that
 * gochem uses, on my own type (which should implement all the relevant gonum interfaces).
 * all the gonum-owned names will start with gn (i.e. RandomFunc becomes gnRandomFunc) so its latter easy to use search and replace to set the
 * correct import path when gonum is implemented (such as gonum.RandomFunc)
 */

const appzero float64 = 0.0000001 //used to correct floating point
//errors. Everything equal or less than this is considered zero.

//The main container, must be able to implement any
//gonum interface.
type CoordMatrix struct {
	*matrix.DenseMatrix
}

//Generate and returns a CoorMatrix from data.
func NewCoords(data []float64, rows, cols int) *CoordMatrix {
	if len(data) < cols*rows {
		panic(NotEnoughElements)
	}
	return &CoordMatrix{matrix.MakeDenseMatrix(data, rows, cols)}

}

//Returns and empty, but not nil, CoordMatrix. It barely allocates memory
func EmptyCoords() *CoordMatrix {
	var a *matrix.DenseMatrix
	return &CoordMatrix{a}

}

//Returns an empty CoordMatrix with the given dimensions
func gnZeros(rows, cols int) *CoordMatrix {
	return &CoordMatrix{matrix.Zeros(rows, cols)}
}

//Returns an identity matrix spanning span cols and rows
func gnEye(span int) *CoordMatrix {
	A := gnZeros(span, span)
	for i := 0; i < span; i++ {
		A.Set(i, i, 1.0)
	}
	return A
}

//This is a facility to sort Eigenvectors/Eigenvalues pairs
//It satisfies the sort.Interface interface.
type eigenpair struct {
	//evecs must have as many rows as evals has elements.
	evecs *CoordMatrix
	evals sort.Float64Slice
}

func (E eigenpair) Less(i, j int) bool {
	return E.evals[i] < E.evals[j]
}
func (E eigenpair) Swap(i, j int) {
	E.evals.Swap(i, j)
	//	E.evecs[i],E.evecs[j]=E.evecs[j],E.evecs[i]
	E.evecs.SwapRows(i, j)
}
func (E eigenpair) Len() int {
	return len(E.evals)
}

//gnEigen wraps the matrix.DenseMatrix.Eigen() function in order to guarantee
//That the eigenvectors and eigenvalues are sorted according to the eigenvalues
//It also guarantees orthonormality and handness. I don't know how many of
//these are already guaranteed by Eig(). Will delete the unneeded parts
//And even this whole function when sure.
func gnEigen(in *CoordMatrix, epsilon float64) (*CoordMatrix, []float64, error) {
	var err error
	if epsilon < 0 {
		epsilon = appzero
	}
	evecsDM, vals, _ := in.Eigen()
	temp := CoordMatrix{evecsDM}
	evecs := &temp
	evals := [3]float64{vals.Get(0, 0), vals.Get(1, 1), vals.Get(2, 2)}
	f := func() { evecs.T(evecs) }
	if err = gnMaybe(gnPanicker(f)); err != nil {
		return nil, nil, err
	}
	eig := eigenpair{evecs, evals[:]}
	sort.Sort(eig)
	//Here I should orthonormalize vectors if needed instead of just complaining.
	//I think orthonormality is guaranteed by  DenseMatrix.Eig() If it is, Ill delete all this
	//If not I'll add ortonormalization routines.
	eigrows, _ := eig.evecs.Dims()
	vectori := EmptyCoords()
	vectorj := EmptyCoords()
	for i := 0; i < eigrows; i++ {
		vectori.RowView(eig.evecs, i)
		for j := i + 1; j < eigrows; j++ {
			vectorj.RowView(eig.evecs, j)
			if math.Abs(vectori.Dot(vectorj)) > epsilon && i != j {
				err = NotOrthogonal //return eig.evecs, evals[:], fmt.Errorf("Vectors not ortogonal!")
			}
		}
		if math.Abs(vectori.Norm(2)-1) > epsilon {
			//Of course I could just normalize the vectors instead of complaining.
			//err= fmt.Errorf("Vectors not normalized %s",err.Error())

		}
	}
	//Checking and fixing the handness of the matrix.This if-else is Jannes idea,
	//I don't really know whether it works.
	eig.evecs.T(eig.evecs)
	if eig.evecs.Det() < 0 {
		eig.evecs.Scale(-1, eig.evecs) //SSC
	} else {
		/*
			eig.evecs.TransposeInPlace()
			eig.evecs.ScaleRow(0,-1)
			eig.evecs.ScaleRow(2,-1)
			eig.evecs.TransposeInPlace()
		*/
		//	fmt.Println("all good, I guess")
	}
	eig.evecs.T(eig.evecs)
	return eig.evecs, eig.evals, err //Returns a slice of evals
}

//Returns the singular value decomposition of matrix A
func gnSVD(A *CoordMatrix) (*CoordMatrix, *CoordMatrix, *CoordMatrix, error) {
	U, s, V, err := A.SVD()
	theU := CoordMatrix{U}
	sigma := CoordMatrix{s}
	theV := CoordMatrix{V}
	return &theU, &sigma, &theV, err

}

//returns a rows,cols matrix filled with gnOnes.
func gnOnes(rows, cols int) *CoordMatrix {
	gnOnes := gnZeros(rows, cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			gnOnes.Set(i, j, 1)
		}
	}
	return gnOnes
}

func gnMul(A, B *CoordMatrix) *CoordMatrix {
	ar, _ := A.Dims()
	_, bc := B.Dims()
	C := gnZeros(ar, bc)
	C.Mul(A, B)
	return C
}

func RowView(A *CoordMatrix, i int) *CoordMatrix {
	B := EmptyCoords()
	B.RowView(A, i)
	return B
}

func gnClone(A *CoordMatrix) *CoordMatrix {
	r, c := A.Dims()
	B := gnZeros(r, c)
	B.Clone(A)
	return B
}

func gnT(A *CoordMatrix) *CoordMatrix {
	r, c := A.Dims()
	B := gnZeros(c, r)
	B.T(A)
	return B
}

//Methods
/* When gonum is ready, all this functions will take a num.Matrix interface as an argument, instead of a
 * CoordMatrix*/

func (F *CoordMatrix) Add(A, B *CoordMatrix) {
	ar, ac := A.Dims()
	br, bc := B.Dims()
	fr, fc := F.Dims()
	if ac != bc || br != ar || ac != fc || ar != fr {
		panic(gnErrShape)
	}
	for i := 0; i < fr; i++ {
		for j := 0; j < fc; j++ {
			F.Set(i, j, A.At(i, j)+B.At(i, j))
		}
	}

}

//AddFloat puts in the receiver a matrix which elements are
//those of matrix A plus the float B.
func (F *CoordMatrix) AddFloat(A *CoordMatrix, B float64) {
	ar, ac := A.Dims()
	if F != A {
		F.Clone(A)
	}
	for i := 0; i < ar; i++ {
		for j := 0; j < ac; j++ {
			F.Set(i, j, A.At(i, j)+B)
		}
	}
}

//AddRow adds the row vector row to each row of the matrix A, putting
//the result on the receiver. Panics if matrices are mismatched.
func (F *CoordMatrix) AddRow(A, row *CoordMatrix) {
	ar, ac := A.Dims()
	rr, rc := row.Dims()
	fr, fc := F.Dims()
	if ac != rc || rr != 1 || ac != fc || ar != fr {
		panic(gnErrShape)
	}
	j := EmptyCoords()
	for i := 0; i < ar; i++ {
		j.RowView(A, i)
		j.Add(j, row)
	}
}

func (F *CoordMatrix) At(A, B int) float64 {
	return F.Get(A, B)
}

func (F *CoordMatrix) Clone(A *CoordMatrix) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if ac != fc || ar != fr {
		panic(gnErrShape)
	}

	for i := 0; i < ar; i++ {
		for j := 0; j < ac; j++ {
			F.Set(i, j, A.At(i, j))
		}

	}

}

//Puts a view of the given row of the matrix on the receiver
func (F *CoordMatrix) ColView(A *CoordMatrix, i int) {
	ar, _ := A.Dims()
	F.View(A, 0, i, ar, 1)
}

func (F *CoordMatrix) DelRow(A *CoordMatrix, i int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if i > ar || fc != ac || fr != (ar-1) {
		panic(gnErrShape)
	}
	tempA1 := EmptyCoords()
	tempF1 := EmptyCoords()
	tempA1.View(A, 0, 0, i, ac)
	tempF1.View(F, 0, 0, i, ac)
	tempF1.Clone(tempA1)
	//now the other part
	tempA2 := EmptyCoords()
	tempF2 := EmptyCoords()
	tempA2.View(A, i+1, 0, ar-i-1, ac) //The magic happens here
	tempF2.View(F, i, 0, ar-i-1, fc)
	tempF2.Clone(tempA2)
}

func (F *CoordMatrix) Dims() (int, int) {
	return F.Rows(), F.Cols()
}

//Dot returns the dot product between 2 vectors or matrices
func (F *CoordMatrix) Dot(B *CoordMatrix) float64 {
	if F.Cols() != B.Cols() || F.Rows() != B.Rows() {
		panic(gnErrShape)
	}
	a, b := F.Dims()
	A := gnZeros(a, b)
	A.MulElem(F, B)
	return A.Sum()
}

//puts the inverse of B in F or panics if F is non-singular.
//its just a dirty minor adaptation from the code in go.matrix from John Asmuth
//it will be replaced by the gonum implementation when the library is ready.
func (F *CoordMatrix) Inv(A *CoordMatrix) {
	//fr,fc:=F.Dims()
	ar, ac := A.Dims()
	if ac != ar {
		panic(gnErrSquare)
	}
	augt, _ := A.Augment(matrix.Eye(ar))
	aug := &CoordMatrix{augt}
	augr, _ := aug.Dims()
	for i := 0; i < augr; i++ {
		j := i
		for k := i; k < augr; k++ {
			if math.Abs(aug.Get(k, i)) > math.Abs(aug.Get(j, i)) {
				j = k
			}
		}
		if j != i {
			aug.SwapRows(i, j)
		}
		if aug.Get(i, i) == 0 {
			panic(gnErrSingular)
		}
		aug.ScaleRow(i, 1.0/aug.Get(i, i))
		for k := 0; k < augr; k++ {
			if k == i {
				continue
			}
			aug.ScaleAddRow(k, i, -aug.Get(k, i))
		}
	}
	F.SubMatrix(aug, 0, ac, ar, ac)
}

//A slightly modified version of John Asmuth's ParalellProduct function.
func (F *CoordMatrix) Mul(A, B *CoordMatrix) {
	if A.Cols() != B.Rows() {
		panic(gnErrShape)
	}
	Arows, Acols := A.Dims()
	_, Bcols := B.Dims()

	if F == nil {
		F = gnZeros(Arows, Bcols) //I don't know if the final API will allow this.
	}

	in := make(chan int)
	quit := make(chan bool)

	dotRowCol := func() {
		for {
			select {
			case i := <-in:
				sums := make([]float64, Bcols)
				for k := 0; k < Acols; k++ {
					for j := 0; j < Bcols; j++ {
						sums[j] += A.At(i, k) * B.At(k, j)
					}
				}
				for j := 0; j < Bcols; j++ {
					F.Set(i, j, sums[j])
				}
			case <-quit:
				return
			}
		}
	}

	threads := 2

	for i := 0; i < threads; i++ {
		go dotRowCol()
	}

	for i := 0; i < Arows; i++ {
		in <- i
	}

	for i := 0; i < threads; i++ {
		quit <- true
	}

	return
}

func (F *CoordMatrix) MulElem(A, B *CoordMatrix) {
	arows, acols := A.Dims()
	brows, bcols := B.Dims()
	frows, fcols := F.Dims()
	if arows != brows || acols != bcols || arows != frows || acols != fcols {
		panic(gnErrShape)
	}
	for i := 0; i < arows; i++ {
		for j := 0; j < acols; j++ {
			F.Set(i, j, A.At(i, j)*B.At(i, j))
		}

	}
}

func (F *CoordMatrix) Norm(i int) float64 {
	//temporary hack
	if i != 2 {
		panic("only 2-norm is implemented")
	}
	return F.TwoNorm()
}

//Puts A**exp on the receiver. This function could probably
//be written in a concurrent way
func (F *CoordMatrix) Pow(A *CoordMatrix, exp float64) {
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

//Returns an array with the data in the ith row of F
func (F *CoordMatrix) Row(i int) []float64 {
	r, c := F.Dims()
	if i >= r {
		panic("Matrix: Requested row out of bounds")
	}
	a := make([]float64, c, c)
	for j := 0; j < c; j++ {
		a[j] = F.At(i, j)
	}
	return a
}

//Puts a view of the given row of the matrix in the receiver
func (F *CoordMatrix) RowView(A *CoordMatrix, i int) {
	_, ac := A.Dims()
	F.View(A, i, 0, 1, ac)
}

//Scale the matrix A by a number i, putting the result in the received.
func (F *CoordMatrix) Scale(i float64, A *CoordMatrix) {
	if A == F { //if A and F points to the same object.
		F.scaleAux(i)
	} else {
		F.Clone(A)
		F.scaleAux(i)
	}
}

func (F *CoordMatrix) scaleAux(factor float64) {
	fr, fc := F.Dims()
	for i := 0; i < fr; i++ {
		for j := 0; j < fc; j++ {
			F.Set(i, j, F.At(i, j)*factor)
		}

	}
}

//ScaleByCol scales each column of matrix A by Col, putting the result
//in the received.
func (F *CoordMatrix) ScaleByCol(A, Col *CoordMatrix) {
	ar, ac := A.Dims()
	cr, cc := Col.Dims()
	fr, fc := F.Dims()
	if ar != cr || cc > 1 || ar != fr || ac != fc {
		panic(gnErrShape)
	}
	if F != A {
		F.Clone(A)
	}
	temp := EmptyCoords()
	for i := 0; i < ac; i++ {
		temp.ColView(F, i)
		temp.MulElem(temp, Col)
	}

}

//ScaleByRow scales each column of matrix A by Col, putting the result
//in the received.
func (F *CoordMatrix) ScaleByRow(A, Row *CoordMatrix) {
	ar, ac := A.Dims()
	rr, rc := Row.Dims()
	fr, fc := F.Dims()
	if ac != rc || rr != 1 || ar != fr || ac != fc {
		panic(gnErrShape)
	}
	if F != A {
		F.Clone(A)
	}
	temp := EmptyCoords()
	for i := 0; i < ac; i++ {
		temp.RowView(F, i)
		temp.MulElem(temp, Row)
	}
}

//When go.matrix is abandoned it is necesary to implement SetMatrix
//SetMatrix()
//Copies A into F aligning A(0,0) with F(i,j)
func (F *CoordMatrix) SetMatrix(i, j int, A *CoordMatrix) {
	fr, fc := F.Dims()
	ar, ac := A.Dims()
	if ar+i > fr || ac+j > fc {
		panic(gnErrShape)
	}
	for l := 0; l < ar; l++ {
		for m := 0; m < ac; m++ {
			F.Set(l+i, m+j, A.Get(l, m))
		}
	}
}

//Returns a matrix contaning all the ith rows of matrix A,
//where i are the numbers in clist. The rows are in the same order
//than the clist.
func (F *CoordMatrix) SomeRows(A *CoordMatrix, clist []int) {
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

//Returns a matrix contaning all the ith rows of matrix A,
//where i are the numbers in clist. The rows are in the same order
//than the clist. Returns an error instead of panicking.
func (F *CoordMatrix) SomeRowsSafe(A *CoordMatrix, clist []int) (err error) {
	f := func() { F.SomeRows(A, clist) }
	return gnMaybe(gnPanicker(f))
}

func (F *CoordMatrix) SetRows(A *CoordMatrix, clist []int) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()
	if ac != fc || fr < len(clist) || ar < len(clist) {
		panic(gnErrShape)
	}
	for key, val := range clist {
		for j := 0; j < ac; j++ {
			F.Set(val, j, A.Get(key, j))
		}
	}
}

//puts in F a matrix consisting in A over B
func (F *CoordMatrix) Stack(A, B *CoordMatrix) {
	Arows, Acols := A.Dims()
	Brows, Bcols := B.Dims()
	Frows, Fcols := F.Dims()

	if Acols != Bcols || Acols != Fcols || Arows+Brows != Frows {
		panic(gnErrShape)
	}

	for i := 0; i < Arows+Brows; i++ {
		for j := 0; j < Acols; j++ {
			if i < Arows {
				F.Set(i, j, A.At(i, j))
			} else {
				F.Set(i, j, B.At(i, j))
			}
		}
	}

	return
}

//Not tested!!!
func (F *CoordMatrix) Sub(A, B *CoordMatrix) {
	B.Scale(-1, B)
	F.Add(A, B)
	B.Scale(-1, B)
}

//not tested
//returns a copy of the submatrix of A starting by the point i,j and
//spanning rows rows and cols columns.
func (F *CoordMatrix) SubMatrix(A *CoordMatrix, i, j, rows, cols int) {
	temp := CoordMatrix{A.GetMatrix(i, j, rows, cols)}
	F.Clone(&temp)
}

//AddRow subtracts the row vector row to each row of the matrix A, putting
//the result on the receiver. Panics if matrices are mismatched.
func (F *CoordMatrix) SubRow(A, row *CoordMatrix) {
	row.Scale(-1, row)
	F.AddRow(A, row)
	row.Scale(-1, row)
}

//Sum returns the sum of all elements in matrix A.
func (F *CoordMatrix) Sum() float64 {
	Rows, Cols := F.Dims()
	var sum float64
	for i := 0; i < Cols; i++ {
		for j := 0; j < Rows; j++ {
			sum += F.Get(j, i)
		}
	}
	return sum
}

//Transpose
func (F *CoordMatrix) T(A *CoordMatrix) {
	ar, ac := A.Dims()
	fr, fc := F.Dims()

	if ar != fc || ac != fr {
		panic(gnErrShape)
	}
	//we do it in a different way if you pass the received as the argument
	//(transpose in place)
	if F == A {
		for i := 0; i < ar; i++ {
			for j := 0; j < i; j++ {
				tmp := A.At(i, j)
				F.Set(i, j, A.At(j, i))
				F.Set(j, i, tmp)
			}

		}
	} else {
		for i := 0; i < ar; i++ {
			for j := 0; j < ac; j++ {
				F.Set(j, i, A.At(i, j))
			}
		}
	}
}

//Unit takes a vector and divides it by its norm
//thus obtaining an unitary vector pointing in the same direction as
//vector.
func (F *CoordMatrix) Unit(A *CoordMatrix) {
	norm := 1.0 / A.Norm(2)
	F.Scale(norm, A)
}

func (F *CoordMatrix) View(A *CoordMatrix, i, j, rows, cols int) {
	*F = CoordMatrix{A.GetMatrix(i, j, rows, cols)}
}

/**These are from the current proposal for gonum, by Dan Kortschak. It will be taken out
 * from here when gonum is implemented. The gn prefix is appended to the names to make them
 * unimported and to allow easy use of search/replace to add the "num" prefix when I change to
 * gonum.**/

// A Panicker is a function that may panic.
type gnPanicker func()

// Maybe will recover a panic with a type matrix.Error from fn, and return this error.
// Any other error is re-panicked.
func gnMaybe(fn gnPanicker) (err error) {
	defer func() {
		if r := recover(); r != nil {
			var ok bool
			if err, ok = r.(gnError); ok {
				return
			}
			panic(r)
		}
	}()
	fn()
	return
}

// Type Error represents matrix package errors. These errors can be recovered by Maybe wrappers.
type gnError string

func (err gnError) Error() string { return string(err) }

const (
	//RM: the first 2 are mine.
	NotOrthogonal        = gnError("matrix: Vectors nor orthogonal")
	NotEnoughElements    = gnError("matrix: not enough elements")
	gnErrIndexOutOfRange = gnError("matrix: index out of range")
	gnErrZeroLength      = gnError("matrix: zero length in matrix definition")
	gnErrRowLength       = gnError("matrix: row length mismatch")
	gnErrColLength       = gnError("matrix: col length mismatch")
	gnErrSquare          = gnError("matrix: expect square matrix")
	gnErrNormOrder       = gnError("matrix: invalid norm order for matrix")
	gnErrSingular        = gnError("matrix: matrix is singular")
	gnErrShape           = gnError("matrix: dimension mismatch")
	gnErrIllegalStride   = gnError("matrix: illegal stride")
	gnErrPivot           = gnError("matrix: malformed pivot list")
)
