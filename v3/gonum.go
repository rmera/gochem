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
	"math"

	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"

	//	"math/cmplx"
	"sort"
)

// Matrix is a set of vectors in 3D space. The underlying implementation varies.
// Within the package it is understood that a "vector" is a row vector, i.e. the
// cartesian coordinates of a point in 3D space. The name of some functions in
// the library reflect this.
type Matrix struct {
	//The main container, must be able to implement most
	//gonum interfaces
	*mat.Dense
}

// Matrix2Dense returns the A Gonum Dense matrix
// associated with A. Changes in one will
// be reflected in t he other
func Matrix2Dense(A *Matrix) *mat.Dense {
	return A.Dense
}

// Dense2Matrix returns a *v3.Matrix
// from a Gonum dense matrix, which has to be
// Nx3
func Dense2Matrix(A *mat.Dense) *Matrix {
	r, c := A.Dims()
	if c != 3 {
		panic(fmt.Sprintf("malformed *mat.Dense matrix to make *v3.Matrix, must be Nx3, is %d x %d", r, c))
	}
	return &Matrix{A}
}

// NewMatrix creates and returns a Matrix with 3 columns from data.
func NewMatrix(data []float64) (*Matrix, error) {
	const cols int = 3
	l := len(data)
	rows := l / cols
	if l%cols != 0 {
		return nil, Error(fmt.Errorf("NewMatrix: Input slice lenght %d not divisible by %d: %d", rows, cols, rows%cols))
	}
	r := mat.NewDense(rows, cols, data)
	return &Matrix{r}, nil
}

// RawSlice returns the underlying []float64 slice for the receiver.
// Changes on either the []float64 or the receiver are expected to
// reflect on the other.
func (F *Matrix) RawSlice() []float64 {
	return F.RawMatrix().Data
}

// Row fills the  dst slice of float64 with the ith row of matrix F and returns it.
// The slice must have the correct size or be nil, in which case a new slice will be created.
// This method is merely a frontend for the mat64.Row function of gonum.
func (F *Matrix) Row(dst []float64, i int) []float64 {
	return mat.Row(dst, i, F.Dense)
}

// Col fills the  dst slice of float64 with the ith col of matrix F and returns it.
// The slice must have the correct size or be nil, in which case a new slice will be created.
// This method is merely a frontend for the mat64.Col function of gonum.
func (F *Matrix) Col(dst []float64, i int) []float64 {
	return mat.Col(dst, i, F.Dense)
}

// ColSlice puts a view of the given col of the matrix on the receiver
func (F *Matrix) ColSlice(i int) *Matrix {
	//	r := new(mat64.Dense)
	Fr, _ := F.Dims()
	r := F.Dense.Slice(0, Fr, i, i+1).(*mat.Dense)
	return &Matrix{r}
}

// ColView a view of the given col of the matrix on the receiver.
// This function is for compatibility with the gonum v1 API
// The older one might be deleted in the future, but, if at all,
// it will take time.
func (F *Matrix) ColView(i int) *Matrix {
	//	r := new(mat64.Dense)
	Fr, _ := F.Dims()
	r := F.Dense.Slice(0, Fr, i, i+1).(*mat.Dense)
	return &Matrix{r}
}

// VecView eturns view of the ith vector of the matrix in the receiver
func (F *Matrix) VecView(i int) *Matrix {
	//r := new(mat64.Dense)
	/*	mr,mc:=F.Caps() /////////////////////////
		j:=0
		k:=1+i
		l:=3
		println("ESTAMOS EN LA BEEE")
		if i < 0{
			println(1)
		}else if mr <= i{
			println(2)
		} else if j < 0{
			println(3)
		} else if mc <= j{
			println(4)
		} else if  k <= i{
			println(5,k,i)
		} else if mr < k {
			println(6)
		} else if l <= j {
			println(7)
		}else if mc < l {
			println(8)
		}
		println("nooo",i,mr,mc)////////////
	*/
	r := F.Dense.Slice(i, i+1, 0, 3).(*mat.Dense)
	return &Matrix{r}
}

// VecSlice slice of the given vector of the matrix in the receiver
// This function is to keep compatibility with the new gonum v1 API
func (F *Matrix) VecSlice(i int) *Matrix {
	//r := new(mat64.Dense)
	r := F.Dense.Slice(i, i+1, 0, 3).(*mat.Dense)
	return &Matrix{r}
}

// View returns a view of F starting from i,j and spanning r rows and
// c columns. Changes in the view are reflected in F and vice-versa
// This view has the wrong signature for the interface mat64.Viewer,
// But the right signatur was not possible to implement. Notice that very little
// memory allocation happens, only a couple of ints and pointers.
func (F *Matrix) View(i, j, r, c int) *Matrix {
	ret := F.Dense.Slice(i, i+r, j, j+c).(*mat.Dense)
	return &Matrix{ret}
}

// Slice returns a view of F starting from i,j and spanning r rows and
// c columns. Changes in the view are reflected in F and vice-versa
// This function is to keep compatibility with the gonum v1 API.
func (F *Matrix) Slice(i, r, j, c int) *Matrix {
	ret := F.Dense.Slice(i, r, j, c).(*mat.Dense)
	return &Matrix{ret}
}

// Sub puts the element-wise subtraction
// of matrices A and B into the receiver
func (F *Matrix) Sub(A, B *Matrix) {
	F.Dense.Sub(A.Dense, B.Dense)
}

// Add puts the element-wise addition
// of matrices A and B into the receiver
func (F *Matrix) Add(A, B *Matrix) {
	F.Dense.Add(A.Dense, B.Dense)
}

// SetMatrix the matrix A in the received starting from the ith row and jth col
// of the receiver.
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

// Mul Wraps mat.Mul to take care of the case when one of the
// arguments is also the received. Since the received is a Matrix,
// the mat64 function could check A (mat64.Dense) vs F (Matrix) and
// it would not know that internally F.Dense==A, hence the need for this function.
func (F *Matrix) Mul(A, B mat.Matrix) {
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

func stupidDot(A, B *Matrix) float64 {
	return A.At(0, 0)*B.At(0, 0) + A.At(0, 1)*B.At(0, 1) + A.At(0, 2)*B.At(0, 2)
}

// Dot gets the dot product between the first row of F and the first row of A. It's a vector dot product,
// to be used with 1-row matrices.
func (F *Matrix) Dot(A *Matrix) float64 {
	//The reason for making Dot ask for a v3.Matrix is that then we can call mat64.Dot with A.Dense, which should make things faster.
	id := mat.NewDense(3, 3, []float64{1, 0, 0, 0, 1, 0, 0, 0, 1}) //Identity matrix
	return mat.Inner(F.Dense.RowView(0), id, A.Dense.RowView(0))
}

// Scale multiplies each element in the matrix A by
// v
func (F *Matrix) Scale(v float64, A *Matrix) {
	F.Dense.Scale(v, A.Dense)
}

// Norm acts as a front-end for the mat64 function.
func (F *Matrix) Norm(i float64) float64 {
	//This used to always return Frobenius norm, no matter what you give as an argument.
	//The argument is there for compatibility (Gonum used to have "0" as the Froebius norm
	//and that was, until recently, still used in goChem. I think I have fixes all these
	//use cases but be mindful of possible bugs.

	return mat.Norm(F.Dense, i)
}

/*
For now I'm not sure we need this wrapper and I do try to keep the to the minimum.

//Sum acts as a wrapper for the mat64 function, for compatibility. It returns the sum of the elements of F.
func (F *Matrix) Sum() float64{
	return mat64.Sum(F.Dense)
}
*/

// Stack puts A stacked over B in the receiver
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

// This is a temporal function. It returns the determinant of a 3x3 matrix. Panics if the matrix is not 3x3.
// It is also defined in the chem package which is not-so-clean.
func det(A mat.Matrix) float64 {
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

// EigenWrap wraps the mat.Eigen structure in order to guarantee
// That the eigenvectors and eigenvalues are sorted according to the eigenvalues
// It also guarantees orthonormality and handness. I don't know how many of
// these are already guaranteed by Eig(). Will delete the unneeded parts
// And even this whole function when sure. The main reason for this function
// Is the compatibiliy with go.matrix. This function should dissapear when we
// have a pure Go blas.
func EigenWrap(in *Matrix, epsilon float64) (*Matrix, []float64, error) {
	if epsilon < 0 {
		epsilon = appzero
	}
	eigen := new(mat.Eigen)
	ok := eigen.Factorize(mat.DenseCopyOf(in.Dense), mat.EigenRight) //Not sure if that DenseCopy is still needed.
	if !ok {
		return nil, nil, Error(fmt.Errorf("error in EigenWrap: mat.Eigen.Factorize"))
	}
	evals_cmp := make([]complex128, 3)  //We only deal with 3-column matrixes in this package
	TempVec := mat.NewCDense(3, 3, nil) /// An allocation. We have to see whether I should try to minimize these. Also, see comment above.
	evals_cmp = eigen.Values(evals_cmp)
	eigen.VectorsTo(TempVec)
	evecsprev := Zeros(3)
	//This is horrible, but, apparently, Gonum just doesn't provide anything to go from CDense to Dense. At least these are guaranteed to be 3x3 matrices
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			scalar := TempVec.At(i, j)
			if imag(scalar) != 0 {
				return nil, nil, Error(fmt.Errorf("mat.EigenFactorize: Found a complex Eigenvector"))
			}
			evecsprev.Set(i, j, real(scalar))
		}
	}
	evals := make([]float64, 3, 3)

	for k, _ := range evals {
		evals[k] = real(evals_cmp[k]) //no check of the thing being real for now.
	}
	evecs := Zeros(3)
	fn := func() { evecs.Copy(evecsprev.T()) }
	err := mat.Maybe(fn)
	if err != nil {
		return nil, nil, Error(fmt.Errorf("EigenWrap: mat.Copy/math.T"))

	}
	//evecs.TCopy(evecs.Dense)
	eig := eigenpair{evecs, evals[:]}
	sort.Sort(eig)
	//Here I should orthonormalize vectors if needed instead of just complaining.
	//I think orthonormality is guaranteed by  DenseMatrix.Eig() If it is, Ill delete all this
	//If not I'll add ortonormalization routines.
	eigrows, _ := eig.evecs.Dims()
	for i := 0; i < eigrows; i++ {
		//vectori := eig.evecs.VecView(i)
		for j := i + 1; j < eigrows; j++ {
			//	vectorj := eig.evecs.VecView(j)
			if math.Abs(eig.evecs.RowView(i).Dot(eig.evecs.RowView(j))) > epsilon && i != j {
				fmt.Println("Dot should be", stupidDot(eig.evecs.RowView(i), eig.evecs.RowView(j)))
				reterr := fmt.Errorf("EigenWrap: Eigenvectors %d and %d nor orthogonal: %v %v. Dot: %5.3f. EigMatrix: %v", i, j, eig.evecs.Dense.RowView(i), eig.evecs.Dense.RowView(j), math.Abs(eig.evecs.RowView(i).Dot(eig.evecs.RowView(j))), eig.evecs)
				return eig.evecs, evals[:], reterr
			}
		}
		if math.Abs(eig.evecs.VecView(1).Norm(2)-1) > epsilon {
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
	}
	//	eig.evecs.TCopy(eig.evecs)
	return eig.evecs, eig.evals, nil //Returns a slice of evals
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

/*
//A wrapper for mat64.Dense.T which returns a Matrix.
func (F *Matrix) Tr () *Matrix{
	Tra:=F.Dense.T()
	if Tra,ok:= Tra.(*mat64.Dense);ok{
		return &Matrix{Tra}
	}else{
		panic("goChem/v3: gonum/matrix/mat64.Dense.T() returned a non mat64.Dense")
	}
}
*/

// I know, premature optimization and so on. It's an internal thing, sue me.
var transposetmp float64

// Tr performs an explicit, in-place tranpose of the receiver.
// it relies in the fact that v3 matrix are all 3D. If the receiver has more than
// 3 rows, the square submatrix of the first 3 rows will be transposed (i.e. no panic or returned error).
// it panics if the receiver has less than 3 rows.
func (F *Matrix) Tr() {
	//This function exists because I can't use the implicit tranpose provided by mat64.Dense.T()
	//which returns a matrix that is not possible to cast into a mat64.Dense
	if F.NVecs() < 3 {
		panic("goChem/v3/Tr: Only 3x3 matrices are allowed for both the argument of Tr(), while the receiver must have 3 rows or more")
	}
	R := F.RawMatrix()
	dataSwitch(R, 0, 1)
	dataSwitch(R, 0, 2)
	dataSwitch(R, 1, 2)
}

// Returns the transpose of F. Does not modify F.
func (F *Matrix) TrRet(toret ...*Matrix) *Matrix {
	if F.NVecs() < 3 {
		panic("goChem/v3/TrRet: The receiver must have 3 rows or more")
	}
	var r *Matrix
	if len(toret) > 0 && toret[0] != nil && toret[0].NVecs() == 3 {
		r = toret[0]
	} else {
		r = Zeros(3)
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r.Set(i, j, F.At(j, i))
		}
	}
	return r
}

// I can only hope this gets inlined
func dataSwitch(R blas64.General, r, c int) {
	transposetmp = R.Data[3*r+c]
	R.Data[3*r+c] = R.Data[3*c+r]
	R.Data[3*c+r] = transposetmp
}

type Error error

// PanicMsg is a message used for panics, even though it does satisfy the error interface.
// for errors use Error.
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
