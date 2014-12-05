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

package chem

import (
	"fmt"
	"github.com/gonum/matrix/mat64"
	"math"
	"sort"
)



//The main container, must be able to implement any
//gonum interface.
//VecMatrix is a set of vectors in 3D space. The underlying implementation varies.
type VecMatrix	struct {
  *mat64.Dense
}

func VecMatrix2Dense(A *VecMatrix) *mat64.Dense {
	return  A.Dense
}

func Dense2VecMatrix(A *mat64.Dense) *VecMatrix {
	return  &VecMatrix{A}
}

//Generate and returns a VecMatrix with 3 columns from data.
func NewVecs(data []float64) (*VecMatrix, error) {
	const cols int = 3
	l := len(data)
	rows := l / cols
	if l%cols != 0 {
		return nil, fmt.Errorf("Input slice lenght %d not divisible by %d: %d", rows, cols, rows%cols)
	}
	r := mat64.NewDense(rows, cols, data)
	return &VecMatrix{r}, nil
}

//Puts a view of the given col of the matrix on the receiver
func (F *VecMatrix) ColView(i int) *VecMatrix {
//	r := new(mat64.Dense)
	Fr, _ := F.Dims()
	r:=F.Dense.View( 0, i, Fr, 1).(*mat64.Dense)
	return &VecMatrix{r}
}

//Returns view of the given vector of the matrix in the receiver
func (F *VecMatrix) VecView(i int) *VecMatrix {
	//r := new(mat64.Dense)
	r:=F.Dense.View( i, 0, 1, 3).(*mat64.Dense)
	return &VecMatrix{r}
}

//View returns a view of F starting from i,j and spanning r rows and
//c columns. Changes in the view are reflected in F and vice-versa
//This view has the wrong signature for the interface mat64.Viewer,
//But the right signatur was not possible to implement. Notice that very little
//memory allocation happens, only a couple of ints and pointers.
func (F *VecMatrix) View(i, j, r, c int) *VecMatrix {
	ret:=F.Dense.View( i, j, r, c).(*mat64.Dense)
	return &VecMatrix{ret}
}

//Puts the matrix A in the received starting from the ith row and jth col
//of the receiver.
func (F *VecMatrix) SetMatrix(i, j int, A *VecMatrix) {
	b := F.RawMatrix()
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

func gnInverse(F *VecMatrix) (*VecMatrix, error) {
	a,err := mat64.Inverse(F.Dense)
	return &VecMatrix{a}, err

}

//Mul Wrapps mat64.Mul to take care of the case when one of the
//argumenst is also the receiver.
/*
func (F *VecMatrix) Mul(A, B mat64.Matrix) {
/*	if F == A {
		A := A.(*VecMatrix)
		F.Dense.Mul(A.Dense, B)
	} else if F == B {
		B := B.(*VecMatrix)
		F.Dense.Mul(A, B.Dense)
	} else {
		F.Dense.Mul(A, B)
	}

	if C, ok := A.(*VecMatrix); ok {
		if D, ok2 := B.(*VecMatrix); ok2 {
				F.Dense.Mul(C.Dense, D.Dense)
			}else{
				F.Dense.Mul(C.Dense,D)
			}

		}else{
			if D, ok2 := B.(*VecMatrix); ok2 {
				F.Dense.Mul(A,D.Dense)
		}else{
			F.Dense.Mul(A,B)
		}
	}
}
*/


//puts A stacked over B in F
func (F *VecMatrix) Stack(A, B *VecMatrix) {
	f := F.RawMatrix()
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


//Some of the following function have an err return type in the signature, but they always return a nil error. This is
//Because of a change in gonum/matrix. The NewDense function used to return error and now panics.
//I do not think it is worth to fix these functions.


//Returns and empty, but not nil, Dense. It barely allocates memory
func emptyDense() (*mat64.Dense, error) {
	a := make([]float64, 0, 0)
	return mat64.NewDense(0, 0, a), nil

}

//Returns an zero-filled Dense with the given dimensions
//It is to be substituted by the Gonum function.
func gnZeros(r, c int) *mat64.Dense {
	f := make([]float64, r*c, r*c)
	return mat64.NewDense(r, c, f)

}

//Returns an identity matrix spanning span cols and rows
func gnEye(span int) *mat64.Dense {
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
				reterr:=VecError(fmt.Sprintln("Eigenvectors ",i,"and",j," not orthogonal. v",i,":",vectori,"\nv",j,":",vectorj,"\nDot:", math.Abs(vectori.Dot(vectorj)),"eigmatrix:", eig.evecs))
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
func (F *VecMatrix) TCopy(A mat64.Matrix) {
	//Somehow the mat64.TCopy method seems to misbehave if I give it a mat64.Matrix.
	//Although I can't see a bug in the mat64.Dense.TCopy function, it seems that if I
	//call it with an A which is not a mat64.Dense, it doesn't work. That is why this wrapper
	//has not been deleted. This seems to be a bug in gochem somehow, not in gonum.
	//to the gonum guys.
	if A, ok := A.(*VecMatrix); ok {
		F.Dense.TCopy(A.Dense)
	} else {
		F.Dense.TCopy(A)
	}
}


//returns a rows,cols matrix filled with gnOnes.
func gnOnes(rows, cols int) *mat64.Dense {
	gnOnes := gnZeros(rows, cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			gnOnes.Set(i, j, 1)
		}
	}
	return gnOnes
}

//The 2 following functions may not even be used.
func gnMul(A, B mat64.Matrix) *mat64.Dense {
	ar, _ := A.Dims()
	_, bc := B.Dims()
	C := gnZeros(ar, bc)
	C.Mul(A, B)
	return C
}

func gnCopy(A mat64.Matrix) *mat64.Dense {
	r, c := A.Dims()
	B := gnZeros(r, c)
	B.Copy(A)
	return B
}

func gnT(A mat64.Matrix) *mat64.Dense {
	r, c := A.Dims()
	B := gnZeros(c, r)
	B.TCopy(A)
	return B
}




//This is a temporal function. It returns the determinant of a 3x3 matrix. Panics if the matrix is not 3x3
func det(A mat64.Matrix) float64 {
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




// A gnPanicker is a function that may panic.
type gnPanicker func()

// Maybe will recover a panic with a type mat64.Error or a VecError from fn, and return this error.
// Any other error is re-panicked. It is a small modification
//Maybe this funciton should be exported.
func gnMaybe(fn gnPanicker) error {
	var err error
	defer func() error {
		if r := recover(); r != nil {
			switch err := r.(type) {
			case VecError:
				return err
			case mat64.Error:
				return err
			default:
			panic(r)
			}
		}
	return nil
	}()
	fn()
	return err
}

// Type Error represents matrix package errors. These errors can be recovered by gnMaybe wrappers.


type VecError string

func (err VecError) Error() string { return string(err) }

const (
	not3xXMatrix      = VecError("goChem/VecMatrix: A VecMatrix should have 3 columns")
	notOrthogonal     = VecError("goChem/VecMatrix: Vectors nor orthogonal")
	notEnoughElements = VecError("goChem/VecMatrix: not enough elements in VecMatrix")
)
