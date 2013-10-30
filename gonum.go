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
	"github.com/gonum/matrix/mat64"
	"github.com/gonum/matrix/mat64/la"
	"math"
	"sort"
)

/*Here I make a -very incomplete- implementation of the gonum api backed by go.matrix, which will enable me to port gochem to gonum.
 * Since the agreement in the gonum community was NOT to build a temporary implementation, I just build the functions that
 * gochem uses, on my own type (which should implement all the relevant gonum interfaces).
 * all the gonum-owned names will start with gn (i.e. RandomFunc becomes gnRandomFunc) so its latter easy to use search and replace to set the
 * correct import path when gonum is implemented (such as gonum.RandomFunc)
 */

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

func ChemDense2VecMatrix(A *ChemDense) *VecMatrix {
	return &VecMatrix{A.Dense}
}

func VecMatrix2ChemDense(A *VecMatrix) *ChemDense {
	return &ChemDense{A.Dense}
}

//Generate and returns a VecMatrix with 3 columns from data.
func NewVecs(data []float64) (*VecMatrix, error) {
	const cols int = 3
	l := len(data)
	rows := l / cols
	//	if l%cols != 0 {
	//		panic(fmt.Sprintf("Input slice lenght %d not divisible by %d: %d", rows, cols, rows%cols))
	//	}
	r, err := mat64.NewDense(rows, cols, data)
	return &VecMatrix{r}, err
}

//Returns a view of the ith Vecinate. Note that the allocation is minimal
func VecView(a *VecMatrix, i int) *VecMatrix {
	ret := a.VecView(i)
	return ret
}


//Puts the matrix A in the received starting from the ith row and jth col
//of the receiver.
func (F *VecMatrix)SetMatrix(i,j int,A *VecMatrix){
	b:=F.BlasMatrix()
	ar,ac:=A.Dims()
	fc:=3
	r:=make([]float64,ac,ac)
	for k:=0;k<ac;k++{
		A.Row(r,k)
		startpoint:=fc*(k+i)+j
		copy(b.Data[startpoint:startpoint+fc],r)
	}
}

func (F *VecMatrix)SwapVecs(i,j int){
	if i>=F.NVecs() || j>=F.NVecs(){
		panic("Indexes out of range")
	}
	rowi := F.Row(nil, i)
	rowj := F.Row(nil, j)
	for k := 0; k < 3; k++ {
		F.Set(i, k, rowj[k])
		F.Set(j, k, rowi[k])
	}
}

//puts A stacked over B in F
func (F *VecMatrix) Stack(A, B *VecMatrix) {
	b := F.BlasMatrix()
	ar, _ := A.Dims()
	br, _ := B.Dims()
	if B.NVecs() < ar+br {
		panic("Not enough space to stack")
	}
	for i := 0; i < ar; i++ {
		A.Row(b.Data[i*3:i*3+3], i)
	}

	for i := ar; i < br; i++ {
		A.Row(b.Data[i*3:i*3+3], i-ar)
	}

}

//func EmptyVecs() *VecMatrix {
//	dens := EmptyDense()
//	return &VecMatrix{dens}
//
//}

//Returns a zero-filled VecMatrix with cos vectors and 3 in the other dimension.
func ZeroVecs(cos int) *VecMatrix {
	const cols int = 3
	dens := gnZeros(cos, cols)
	return &VecMatrix{dens.Dense}
}

//Just a dense matrix to allow different implementations of gonum.
//This might change to Dense at some point, but the API change should be barely noticeable.
type ChemDense struct {
	*mat64.Dense
}

func NewChemDense(data []float64, r, c int) (*ChemDense, error) {
	d, err := mat64.NewDense(r, c, data)
	return &ChemDense{d}, err
}

//Returns and empty, but not nil, Dense. It barely allocates memory
func emptyDense() (*mat64.Dense, error) {
	a := make([]float64, 0, 0)
	return mat64.NewDense(0, 0, a)

}

//Returns an zero-filled Dense with the given dimensions
//It is to be substituted by the Gonum function.
func gnZeros(r, c int) *ChemDense {
	f := make([]float64, r*c, r*c)
	ret, _ := mat64.NewDense(r, c, f)
	return &ChemDense{ret}

}

//Returns an identity matrix spanning span cols and rows
func gnEye(span int) *ChemDense {
	A := gnZeros(span, span)
	for i := 0; i < span; i++ {
		A.Set(i, i, 1.0)
	}
	return A
}

//func Eye(span int) *ChemDense {
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
	E.evecs.SwapVecs(i,j)

}
func (E eigenpair) Len() int {
	return len(E.evals)
}

//gnEigen wraps the matrix.DenseMatrix.Eigen() function in order to guarantee
//That the eigenvectors and eigenvalues are sorted according to the eigenvalues
//It also guarantees orthonormality and handness. I don't know how many of
//these are already guaranteed by Eig(). Will delete the unneeded parts
//And even this whole function when sure.
func gnEigen(in *VecMatrix, epsilon float64) (*VecMatrix, []float64, error) {
	var err error
	if epsilon < 0 {
		epsilon = appzero
	}
	evals, _, vecs := la.Eigen(in.Dense, epsilon)
	evecs := &VecMatrix{vecs}
	//evals := [3]float64{vals.At(0, 0), vals.At(1, 1), vals.At(2, 2)} //go.matrix specific code here.
	f := func() { evecs.TCopy(evecs) }
	if err = gnMaybe(gnPanicker(f)); err != nil {
		return nil, nil, err
	}
	eig := eigenpair{evecs, evals[:]}
	sort.Sort(eig)
	//Here I should orthonormalize vectors if needed instead of just complaining.
	//I think orthonormality is guaranteed by  DenseMatrix.Eig() If it is, Ill delete all this
	//If not I'll add ortonormalization routines.
	eigrows, _ := eig.evecs.Dims()
	for i := 0; i < eigrows; i++ {
		vectori := eig.evecs.ColView(i)
		for j := i + 1; j < eigrows; j++ {
			vectorj := eig.evecs.ColView(j)
			if math.Abs(vectori.Dot(vectorj)) > epsilon && i != j {
				return eig.evecs, evals[:], NotOrthogonal
			}
		}
		if math.Abs(vectori.Norm(0)-1) > epsilon {
			//Of course I could just normalize the vectors instead of complaining.
			//err= fmt.Errorf("Vectors not normalized %s",err.Error())

		}
	}
	//Checking and fixing the handness of the matrix.This if-else is Jannes idea,
	//I don't really know whether it works.
	eig.evecs.TCopy(eig.evecs)
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
	eig.evecs.TCopy(eig.evecs)
	return eig.evecs, eig.evals, err //Returns a slice of evals
}

//Returns the singular value decomposition of matrix A
func gnSVD(A *ChemDense) (*ChemDense, *ChemDense, *ChemDense) {
	s,U, V := la.SVD(A.Dense, appzero, appzero, true, true) //I am not sure that the second appzero is appropiate
	//make sigma a matrix
	lens:=len(s)
	Sigma, _ := mat64.NewDense(lens, lens, make([]float64, lens*lens)) //the slice is hardcoded, no error
	for i := 0; i < lens; i++ {
		Sigma.Set(i, i, s[i])
	}
	return &ChemDense{U}, &ChemDense{Sigma}, &ChemDense{V}

}

//returns a rows,cols matrix filled with gnOnes.
func gnOnes(rows, cols int) *ChemDense {
	gnOnes := gnZeros(rows, cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			gnOnes.Set(i, j, 1)
		}
	}
	return gnOnes
}

func gnMul(A, B mat64.Matrix) *ChemDense {
	ar, _ := A.Dims()
	_, bc := B.Dims()
	C := gnZeros(ar, bc)
	C.Mul(A, B)
	return C
}

func gnClone(A mat64.Matrix) *ChemDense {
	r, c := A.Dims()
	B := gnZeros(r, c)
	B.Clone(A)
	return B
}

func gnT(A mat64.Matrix) *ChemDense {
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
	Not3xXMatrix      = gnError("matrix: The other dimmension should be 3")
	NotOrthogonal     = gnError("matrix: Vectors nor orthogonal")
	NotEnoughElements = gnError("matrix: not enough elements")
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
