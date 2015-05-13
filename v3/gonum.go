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

//gonum.go contains most of what is needed for handling the gonum/mat64 types and facilities.
//At this point the name is mostly historical: It used to be the only file importing gonum,
//Now that gomatrix support is discontinued, there is a tighter integration with gonum and
//other files import mat64.

//All the *Vec functions will operate/produce column or row vectors depending on whether the matrix underlying Dense
//is row or column major.

package v3

import (
	"fmt"
	"github.com/gonum/matrix/mat64"
	"math"
	"sort"
)

//The main container, must be able to implement any
//gonum interface.
//Matrix is a set of vectors in 3D space. The underlying implementation varies.
//Within the package it is understood that a "vector" is a row vector, i.e. the
//cartesian coordinates of a point in 3D space. The name of some funcitions in
//the library reflect this.
type Matrix struct {
	*mat64.Dense
}

func Matrix2Dense(A *Matrix) *mat64.Dense {
	return A.Dense
}

func Dense2Matrix(A *mat64.Dense) *Matrix {
	return &Matrix{A}
}

//Generate and returns a Matrix with 3 columns from data.
func NewMatrix(data []float64) (*Matrix, error) {
	const cols int = 3
	l := len(data)
	rows := l / cols
	if l%cols != 0 {
		return nil, Error{fmt.Sprintf("Input slice lenght %d not divisible by %d: %d", rows, cols, rows%cols), []string{"NewMatrix"}, true}
	}
	r := mat64.NewDense(rows, cols, data)
	return &Matrix{r}, nil
}

//Puts a view of the given col of the matrix on the receiver
func (F *Matrix) ColView(i int) *Matrix {
	//	r := new(mat64.Dense)
	Fr, _ := F.Dims()
	r := F.Dense.View(0, i, Fr, 1).(*mat64.Dense)
	return &Matrix{r}
}

//Returns view of the given vector of the matrix in the receiver
func (F *Matrix) VecView(i int) *Matrix {
	//r := new(mat64.Dense)
	r := F.Dense.View(i, 0, 1, 3).(*mat64.Dense)
	return &Matrix{r}
}

//View returns a view of F starting from i,j and spanning r rows and
//c columns. Changes in the view are reflected in F and vice-versa
//This view has the wrong signature for the interface mat64.Viewer,
//But the right signatur was not possible to implement. Notice that very little
//memory allocation happens, only a couple of ints and pointers.
func (F *Matrix) View(i, j, r, c int) *Matrix {
	ret := F.Dense.View(i, j, r, c).(*mat64.Dense)
	return &Matrix{ret}
}

//Puts the matrix A in the received starting from the ith row and jth col
//of the receiver.
func (F *Matrix) SetMatrix(i, j int, A *Matrix) {
	b := F.RawMatrix()
	ar, ac := A.Dims()
	fc := 3
	if ar+i > F.NVecs() || ac+j > fc {
		panic(ErrShape)
	}
	r := make([]float64, ac, ac)
	for k := 0; k < ar; k++ {
		A.Row(r, k)
		startpoint := fc*(k+i) + j
		copy(b.Data[startpoint:startpoint+fc], r)
	}
}

//Mul Wrapps mat64.Mul to take care of the case when one of the
//arguments is also the received. Since the received is a Matrix,
//the mat64 function could check A (mat64.Dense) vs F (Matrix) and
//it would not know that internally F.Dense==A, hence the need for this function.
func (F *Matrix) Mul(A, B mat64.Matrix) {
	if F == A {
		A := A.(*Matrix)
		F.Dense.Mul(A.Dense, B)
	} else if F == B {
		B := B.(*Matrix)
		F.Dense.Mul(A, B.Dense)
	} else {
		F.Dense.Mul(A, B)
	}

	/*
		if C, ok := A.(*Matrix); ok {
			if D, ok2 := B.(*Matrix); ok2 {
					F.Dense.Mul(C.Dense, D.Dense)
				}else{
					F.Dense.Mul(C.Dense,D)
				}

			}else{
				if D, ok2 := B.(*Matrix); ok2 {
					F.Dense.Mul(A,D.Dense)
			}else{
				F.Dense.Mul(A,B)
			}
		}
	*/
}

//puts A stacked over B in F
func (F *Matrix) Stack(A, B *Matrix) {
	f := F.RawMatrix()
	ar, _ := A.Dims()
	br, _ := B.Dims()
	if F.NVecs() < ar+br {
		panic(ErrShape)
	}
	for i := 0; i < ar; i++ {
		A.Row(f.Data[i*3:i*3+3], i)
	}

	for i := ar; i < ar+br; i++ {
		B.Row(f.Data[i*3:i*3+3], i-ar)
	}

}

//func EmptyVecs() *Matrix {
//	dens := EmptyDense()
//	return &Matrix{dens}
//
//}

//Some of the following function have an err return type in the signature, but they always return a nil error. This is
//Because of a change in gonum/matrix. The NewDense function used to return error and now panics.
//I do not think it is worth to fix these functions.

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

//This is a temporal function. It returns the determinant of a 3x3 matrix. Panics if the matrix is not 3x3.
//It is also defined in the chem package which is not-so-clean.
func det(A mat64.Matrix) float64 {
	r, c := A.Dims()
	if r != 3 || c != 3 {
		panic(ErrDeterminant)
	}
	return (A.At(0, 0)*(A.At(1, 1)*A.At(2, 2)-A.At(2, 1)*A.At(1, 2)) - A.At(1, 0)*(A.At(0, 1)*A.At(2, 2)-A.At(2, 1)*A.At(0, 2)) + A.At(2, 0)*(A.At(0, 1)*A.At(1, 2)-A.At(1, 1)*A.At(0, 2)))
}

type eigenpair struct {
	//evecs must have as many rows as evals has elements.
	evecs *Matrix
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
func EigenWrap(in *Matrix, epsilon float64) (*Matrix, []float64, error) {
	err := Error{string(ErrEigen), []string{"EigenWrap"}, true}

	if epsilon < 0 {
		epsilon = appzero
	}
	efacs := mat64.Eigen(mat64.DenseCopyOf(in.Dense), epsilon)
	evecs := &Matrix{efacs.V}
	evalsmat := efacs.D()
	d, _ := evalsmat.Dims()
	evals := make([]float64, d, d)
	for k, _ := range evals {
		evals[k] = evalsmat.At(k, k)
	}
	//	fmt.Println("evecs fresh", evecs) ///////
	//evals := [3]float64{vals.At(0, 0), vals.At(1, 1), vals.At(2, 2)} //go.matrix specific code here.
	f := func() { evecs.TCopy(evecs) }
	if err2 := mat64.Maybe(mat64.Panicker(f)); err2 != nil {
		return nil, nil, Error{err.Error(), []string{"EigenWrap"}, true}
	}
	//evecs.TCopy(evecs.Dense)
	//	fmt.Println("evecs presort", evecs) /////////
	eig := eigenpair{evecs, evals[:]}
	//fmt.Println("EVEECS", evecs) /////////////////////////
	sort.Sort(eig)
	//Here I should orthonormalize vectors if needed instead of just complaining.
	//I think orthonormality is guaranteed by  DenseMatrix.Eig() If it is, Ill delete all this
	//If not I'll add ortonormalization routines.
	eigrows, _ := eig.evecs.Dims()
	//	fmt.Println("evecs", eig.evecs) /////////
	for i := 0; i < eigrows; i++ {
		vectori := eig.evecs.VecView(i)
		for j := i + 1; j < eigrows; j++ {
			vectorj := eig.evecs.VecView(j)
			if math.Abs(vectori.Dot(vectorj)) > epsilon && i != j {
				reterr := Error{fmt.Sprintln("Eigenvectors ", i, "and", j, " not orthogonal. v", i, ":", vectori, "\nv", j, ":", vectorj, "\nDot:", math.Abs(vectori.Dot(vectorj)), "eigmatrix:", eig.evecs), []string{"EigenWrap"}, true}
				return eig.evecs, evals[:], reterr
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

/*
//Returns the singular value decomposition of matrix A
func gnSVD(A *mat64.Dense) ( *mat64.Dense,*mat64.Dense,*mat64.Dense) {
	facts := mat64.SVD(A.Dense, appzero, math.SmallestNonzeroFloat64, true, true) //I am not sure that the second appzero is appropiate
	//make sigma a matrix
	//	lens:=len(s)
	//	Sigma, _ := mat64.NewDense(lens, lens, make([]float64, lens*lens)) //the slice is hardcoded, no error
	//	for i := 0; i < lens; i++ {
	//		Sigma.Set(i, i, s[i])
	//	}
	return &chemDense{facts.U}, &chemDense{facts.S()}, &chemDense{facts.V}

}
*/

//Just a wrapper for the mat64.Dense.TCopy method
func (F *Matrix) TCopy(A mat64.Matrix) {
	//Somehow the mat64.TCopy method seems to misbehave if I give it a mat64.Matrix.
	//Although I can't see a bug in the mat64.Dense.TCopy function, it seems that if I
	//call it with an A which is not a mat64.Dense, it doesn't work. That is why this wrapper
	//has not been deleted. This seems to be a bug in gochem somehow, not in gonum.
	if A, ok := A.(*Matrix); ok {
		F.Dense.TCopy(A.Dense)
	} else {
		F.Dense.TCopy(A)
	}
}

//Errors

//the same as chem.Error but avoid circular import.
type errorInt interface {
	Error() string
	Critical() bool
	Decorate(string) []string
}

type Error struct {
	message  string
	deco     []string
	critical bool
}

//Error returns a string with an error message.
func (err Error) Error() string {
	return fmt.Sprintf("%s", err.message)
}

//Decorate will add the dec string to the decoration slice of strings of the error,
//and return the resulting slice.
func (err Error) Decorate(dec string) []string {
	err.deco = append(err.deco, dec)
	return err.deco
}

//Critical return whether the error is critical or it can be ifnored
func (err Error) Critical() bool { return err.critical }

//errDecorate is a helper function that asserts that the error is
//implements chem.Error and decorates the error with the caller's name before returning it.
//if used with a non-chem.Error error, it will cause a panic.
func errDecorate(err error, caller string) error {
	err2 := err.(errorInt) //I know that is the type returned byt initRead
	err2.Decorate(caller)
	return err2
}

//PanicMsg is a message used for panics, even though it does satisfy the error interface.
//for errors use Error.
type PanicMsg string

func (v PanicMsg) Error() string { return string(v) }

//Here I use a few of the gonum/mat64 messages for compatibility (Copyright (c) The gonum authors). I assumed here that this is small enough not to require
//messing about with licenses, but of course I don't intend any copyright impringement and I will set somethign up if contacted.

const (
	ErrNotXx3Matrix      = PanicMsg("goChem/v3: A VecMatrix should have 3 columns")
	ErrNoCrossProduct    = PanicMsg("goChem/v3: Invalid matrix for cross product")
	ErrNotOrthogonal     = PanicMsg("goChem/v3: Vectors nor orthogonal")
	ErrNotEnoughElements = PanicMsg("goChem/v3: not enough elements in Matrix")
	ErrGonum             = PanicMsg("goChem/v3: Error in gonum function")
	ErrEigen             = PanicMsg("goChem/v3: Can't obtain eigenvectors/eigenvalues of given matrix")
	ErrDeterminant       = PanicMsg("goChem/v3: Determinants are only available for 3x3 matrices")
	ErrShape             = PanicMsg("goChem/v3: Dimension mismatch")
	ErrIndexOutOfRange   = PanicMsg("mat64: index out of range")
)
