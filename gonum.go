// +build !purego

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

//All the *Vec functions will operate/produce column or row vectors depending on whether the matrix underlying Dense
//is row or column major.

package chem

import (
	"fmt"
	"github.com/gonum/blas/cblas"
	"github.com/gonum/matrix/mat64"
	"math"
	"sort"
)

//For now this is here as we do not have other blas engine options.
//When we do, there will be several files with different inits,
//That will be chosen with compiler flags.
func init() {
	mat64.Register(cblas.Blas{})
}

//INTERFACES  This part is from GONUM, copyright, the gonum authors.

//The main container, must be able to implement any
//gonum interface.
//VecMatrix is a set of vectors in 3D space. The underlying implementation varies.
type VecMatrix struct {
	*mat64.Dense
}

func VecMatrix2Dense(A *VecMatrix) *mat64.Dense {
	return A.Dense
}

func Dense2VecMatrix(A *mat64.Dense) *VecMatrix {
	return &VecMatrix{A}
}

//Generate and returns a VecMatrix with 3 columns from data.
func NewVecs(data []float64) (*VecMatrix, error) {
	const cols int = 3
	l := len(data)
	rows := l / cols
	if l%cols != 0 {
		return nil, fmt.Errorf("Input slice lenght %d not divisible by %d: %d", rows, cols, rows%cols)
	}
	r, err := mat64.NewDense(rows, cols, data)
	return &VecMatrix{r}, err
}

//Puts a view of the given col of the matrix on the receiver
func (F *VecMatrix) ColView(i int) *VecMatrix {
	r := new(mat64.Dense)
	*r = *F.Dense
	Fr, _ := F.Dims()
	r.View(0, i, Fr, 1)
	return &VecMatrix{r}
}

//Returns view of the given vector of the matrix in the receiver
func (F *VecMatrix) VecView(i int) *VecMatrix {
	r := new(mat64.Dense)
	*r = *F.Dense
	r.View(i, 0, 1, 3)
	return &VecMatrix{r}
}

//View returns a view of F starting from i,j and spanning r rows and
//c columns. Changes in the view are reflected in F and vice-versa
//This view has the wrong signature for the interface mat64.Viewer,
//But the right signatur was not possible to implement. Notice that very little
//memory allocation happens, only a couple of ints and pointers.
func (F *VecMatrix) View(i, j, r, c int) *VecMatrix {
	ret := new(mat64.Dense)
	*ret = *F.Dense
	ret.View(i, j, r, c)
	return &VecMatrix{ret}
}

//Puts the matrix A in the received starting from the ith row and jth col
//of the receiver.
func (F *VecMatrix) SetMatrix(i, j int, A *VecMatrix) {
	b := F.BlasMatrix()
	ar, ac := A.Dims()
	fc := 3
	if ar+i > F.NVecs() || ac+j > fc {
		panic("SetMatrix: Dimmension mismatch")
	}
	r := make([]float64, ac, ac)
	for k := 0; k < ar; k++ {
		A.Row(r, k)
		startpoint := fc*(k+i) + j
		copy(b.Data[startpoint:startpoint+fc], r)
	}
}

//The following functions are here because directly calling the method
//on the embeded mat64.Dense gives wrong behavior.

/*

func (F *VecMatrix) Copy(A Matrix) {
	switch A := A.(type) {
	case *VecMatrix:
		F.Dense.Copy(A.Dense)
	case *chemDense:
		F.Dense.Copy(A.Dense)
	default:
		F.Dense.Copy(A)
	}
}


func (F *VecMatrix)Add(A *VecMatrix, B Matrix){
	switch B := B.(type) {
	case *VecMatrix:
		F.Dense.Add(A.Dense,B.Dense)
	case *chemDense:
		F.Dense.Add(A.Dense,B.Dense)
	default:
		F.Dense.Add(A.Dense,B)
	}
}

func (F *VecMatrix)Scale(i float64, A Matrix){
	switch A := A.(type) {
	case *VecMatrix:
		F.Dense.Scale(i,A.Dense)
	case *chemDense:
		F.Dense.Scale(i,A.Dense)
	default:
		F.Dense.Scale(i,A)
	}
}


func (F *VecMatrix) Sub(A *VecMatrix, B Matrix){
	switch B := B.(type) {
	case *VecMatrix:
		F.Dense.Sub(A.Dense,B.Dense)
	case *chemDense:
		F.Dense.Sub(A.Dense,B.Dense)
	default:
		F.Dense.Sub(A.Dense,B)
	}
}


*/

func gnInverse(F *VecMatrix) *VecMatrix {
	a := mat64.Inverse(F.Dense)
	return &VecMatrix{a}

}

//Mul Wrapps mat64.Mul to take care of the case when one of the
//argumenst is also the receiver.
func (F *VecMatrix) Mul(A, B Matrix) {
	if F == A {
		A := A.(*VecMatrix)
		F.Dense.Mul(A.Dense, B)
	} else if F == B {
		B := B.(*VecMatrix)
		F.Dense.Mul(A, B.Dense)
	} else {
		F.Dense.Mul(A, B)
	}
	/*
		if C, ok := A.(*VecMatrix); ok {
			switch B := B.(type) {
			case *VecMatrix:
				F.Dense.Mul(C.Dense, B.Dense)
			case *chemDense:
				F.Dense.Mul(C.Dense, B.Dense)
			default:
				F.Dense.Mul(C.Dense, B)
			}
		} else if C,ok:=A.(*chemDense); ok {
			switch B := B.(type) {
			case *VecMatrix:
				F.Dense.Mul(C.Dense, B.Dense)
			case *chemDense:
				F.Dense.Mul(C.Dense, B.Dense)
			default:
				F.Dense.Mul(C.Dense, B)
			}
		} else {
			F.Dense.Mul(A, B)
		}
	*/
}

//puts A stacked over B in F
func (F *VecMatrix) Stack(A, B *VecMatrix) {
	f := F.BlasMatrix()
	ar, _ := A.Dims()
	br, _ := B.Dims()
	if F.NVecs() < ar+br {
		panic("Not enough space to stack")
	}
	for i := 0; i < ar; i++ {
		A.Row(f.Data[i*3:i*3+3], i)
	}

	for i := ar; i < ar+br; i++ {
		B.Row(f.Data[i*3:i*3+3], i-ar)
	}

}

//func EmptyVecs() *VecMatrix {
//	dens := EmptyDense()
//	return &VecMatrix{dens}
//
//}

//Just a dense matrix to allow different implementations of gonum.
//This might change to Dense at some point, but the API change should be barely noticeable.
type chemDense struct {
	*mat64.Dense
}

func newchemDense(data []float64, r, c int) (*chemDense, error) {
	d, err := mat64.NewDense(r, c, data)
	return &chemDense{d}, err
}

//Returns and empty, but not nil, Dense. It barely allocates memory
func emptyDense() (*mat64.Dense, error) {
	a := make([]float64, 0, 0)
	return mat64.NewDense(0, 0, a)

}

//Returns an zero-filled Dense with the given dimensions
//It is to be substituted by the Gonum function.
func gnZeros(r, c int) *chemDense {
	f := make([]float64, r*c, r*c)
	ret, _ := mat64.NewDense(r, c, f)
	return &chemDense{ret}

}

//Returns an identity matrix spanning span cols and rows
func gnEye(span int) *chemDense {
	A := gnZeros(span, span)
	for i := 0; i < span; i++ {
		A.Set(i, i, 1.0)
	}
	return A
}

//func Eye(span int) *chemDense {
//	return gnEye(span)
//}

//Some temporary support function.
//func Eigen(in *Dense, epsilon float64) (*Dense, []float64, error) {
//	i, j, k := gnEigen(in, epsilon)
//	return i, j, k
//}

//This is a facility to sort Eigenvectors/Eigenvalues pairs
//It satisfies the sort.Interface interface.
type eigenpair struct {
	//evecs must have as many rows as evals has elements.
	evecs *VecMatrix
	evals sort.Float64Slice
}

func (E eigenpair) Less(i, j int) bool {
	return E.evals[i] < E.evals[j]
}
func (E eigenpair) Swap(i, j int) {
	E.evals.Swap(i, j)
	//	E.evecs[i],E.evecs[j]=E.evecs[j],E.evecs[i]
	E.evecs.SwapVecs(i, j)

}
func (E eigenpair) Len() int {
	return len(E.evals)
}

//gnEigen wraps the matrix.DenseMatrix.Eigen() function in order to guarantee
//That the eigenvectors and eigenvalues are sorted according to the eigenvalues
//It also guarantees orthonormality and handness. I don't know how many of
//these are already guaranteed by Eig(). Will delete the unneeded parts
//And even this whole function when sure. The main reason for this function
//Is the compatibiliy with go.matrix. This function should dissapear when we
//have a pure Go blas.
func EigenWrap(in *VecMatrix, epsilon float64) (*VecMatrix, []float64, error) {
	var err error
	if epsilon < 0 {
		epsilon = appzero
	}
	efacs := mat64.Eigen(mat64.DenseCopyOf(in.Dense), epsilon)
	evecs := &VecMatrix{efacs.V}
	evalsmat := efacs.D()
	d, _ := evalsmat.Dims()
	evals := make([]float64, d, d)
	for k, _ := range evals {
		evals[k] = evalsmat.At(k, k)
	}
	//	fmt.Println("evecs fresh", evecs) ///////
	//evals := [3]float64{vals.At(0, 0), vals.At(1, 1), vals.At(2, 2)} //go.matrix specific code here.
	f := func() { evecs.TCopy(evecs) }
	if err = mat64.Maybe(mat64.Panicker(f)); err != nil {
		return nil, nil, err
	}
	//evecs.TCopy(evecs.Dense)
	//	fmt.Println("evecs presort", evecs) /////////
	eig := eigenpair{evecs, evals[:]}

	sort.Sort(eig)
	//Here I should orthonormalize vectors if needed instead of just complaining.
	//I think orthonormality is guaranteed by  DenseMatrix.Eig() If it is, Ill delete all this
	//If not I'll add ortonormalization routines.
	eigrows, _ := eig.evecs.Dims()
	//	fmt.Println("evecs", eig.evecs) /////////
	for i := 0; i < eigrows; i++ {
		vectori := eig.evecs.RowView(i)
		for j := i + 1; j < eigrows; j++ {
			vectorj := eig.evecs.RowView(j)
			if math.Abs(vectori.Dot(vectorj)) > epsilon && i != j {
				fmt.Println("FAAAAILL", eig.evecs, i, j, math.Abs(vectori.Dot(vectorj)), vectori, vectorj)
				return eig.evecs, evals[:], notOrthogonal
			}
		}
		if math.Abs(vectori.Norm(0)-1) > epsilon {
			//Of course I could just normalize the vectors instead of complaining.
			//err= fmt.Errorf("Vectors not normalized %s",err.Error())

		}
	}
	//Checking and fixing the handness of the matrix.This if-else is Jannes idea,
	//I don't really know whether it works.
	//	eig.evecs.TCopy(eig.evecs)
	if det(eig.evecs) < 0 { //Right now, this will fail if the matrix is not 3x3 (30/10/2013)
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
	//	eig.evecs.TCopy(eig.evecs)
	return eig.evecs, eig.evals, err //Returns a slice of evals
}

//Returns the singular value decomposition of matrix A
func gnSVD(A *chemDense) (*chemDense, *chemDense, *chemDense) {
	facts := mat64.SVD(A.Dense, appzero, math.SmallestNonzeroFloat64, true, true) //I am not sure that the second appzero is appropiate
	//make sigma a matrix
	//	lens:=len(s)
	//	Sigma, _ := mat64.NewDense(lens, lens, make([]float64, lens*lens)) //the slice is hardcoded, no error
	//	for i := 0; i < lens; i++ {
	//		Sigma.Set(i, i, s[i])
	//	}
	return &chemDense{facts.U}, &chemDense{facts.S()}, &chemDense{facts.V}

}

func (F *VecMatrix) TCopy(A Matrix) {
	if A, ok := A.(*VecMatrix); ok {
		F.Dense.TCopy(A.Dense)
	} else {
		F.Dense.TCopy(A)
	}
}

//returns a rows,cols matrix filled with gnOnes.
func gnOnes(rows, cols int) *chemDense {
	gnOnes := gnZeros(rows, cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			gnOnes.Set(i, j, 1)
		}
	}
	return gnOnes
}

func gnMul(A, B mat64.Matrix) *chemDense {
	ar, _ := A.Dims()
	_, bc := B.Dims()
	C := gnZeros(ar, bc)
	C.Mul(A, B)
	return C
}

func gnCopy(A mat64.Matrix) *chemDense {
	r, c := A.Dims()
	B := gnZeros(r, c)
	B.Copy(A)
	return B
}

func gnT(A mat64.Matrix) *chemDense {
	r, c := A.Dims()
	B := gnZeros(c, r)
	B.TCopy(A)
	return B
}

//This is a temporal function. It returns the determinant of a 3x3 matrix. Panics if the matrix is not 3x3
func det(A Matrix) float64 {
	r, c := A.Dims()
	if r != 3 || c != 3 {
		panic("Determinants are for now only available for 3x3 matrices")
	}
	return (A.At(0, 0)*(A.At(1, 1)*A.At(2, 2)-A.At(2, 1)*A.At(1, 2)) - A.At(1, 0)*(A.At(0, 1)*A.At(2, 2)-A.At(2, 1)*A.At(0, 2)) + A.At(2, 0)*(A.At(0, 1)*A.At(1, 2)-A.At(1, 1)*A.At(0, 2)))
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
	//RM
	not3xXMatrix      = gnError("matrix: The other dimmension should be 3")
	notOrthogonal     = gnError("matrix: Vectors nor orthogonal")
	notEnoughElements = gnError("matrix: not enough elements")
	//end RM
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
