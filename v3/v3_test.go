/*
 * coord_test.go
 *
 * Copyright 2013 Raul Mera <rmera@zinc>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

package v3

import (
	"fmt"
	"testing"

	"gonum.org/v1/gonum/mat"
)

// Returns an zero-filled Dense with the given dimensions
// It is to be substituted by the Gonum function.
func gnZeros(r, c int) *mat.Dense {
	f := make([]float64, r*c, r*c)
	return mat.NewDense(r, c, f)

}

// Returns an identity matrix spanning span cols and rows
func gnEye(span int) *mat.Dense {
	A := gnZeros(span, span)
	for i := 0; i < span; i++ {
		A.Set(i, i, 1.0)
	}
	return A
}

func TestGeo(Te *testing.T) {
	a := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9}
	A, err := NewMatrix(a)
	if err != nil {
		Te.Error(err)
	}
	ar, ac := A.Dims()
	T := Zeros(ar)
	T.Copy(A)
	T.T()
	B := gnEye(ar)
	//B.Copy(A)
	T.Mul(A, B)
	E := Zeros(ar)
	E.MulElem(A, B)
	fmt.Println(T, "\n", T, "\n", A, "\n", B, "\n", ar, ac)
	View := A.VecView(1)
	View.Set(0, 0, 100)
	fmt.Println("View\n", A, "\n", View)

}

func TestSomeVecs(Te *testing.T) {
	a := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}
	A, err := NewMatrix(a)
	if err != nil {
		Te.Error(err)
	}
	B := Zeros(3) //We should cause an error by modifying this.
	cind := []int{1, 3, 5}
	err = B.SomeVecsSafe(A, cind)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println(A, "\n", B)
	B.Set(1, 1, 55)
	B.Set(2, 2, 66)
	fmt.Println("Changed B")
	fmt.Println(A, "\n", B)
	A.SetVecs(B, cind)
	fmt.Println("Now A should see changes in B")
	fmt.Println(A, "\n", B)
}

func TestScale(Te *testing.T) {
	a := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}
	A, err := NewMatrix(a)
	if err != nil {
		Te.Error(err)
	}
	B := Zeros(6)
	A.Scale(3, A)
	B.Scale(2, A)
	fmt.Println(A, "\n", B)
	Row, err := NewMatrix([]float64{10, 20, 30})
	if err != nil {
		Te.Error(err)
	}
	A.AddVec(A, Row)
	fmt.Println("Additions")
	fmt.Println(A)
	A.SubVec(A, Row)
	fmt.Println(A, A.NVecs(), B.NVecs())
	//	B.Pow(A, 2)
	//	fmt.Println("Squared", A, "\n", B)
	b := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9}
	S, err := NewMatrix(b)
	if err != nil {
		Te.Error(err)
	}
	row, err := NewMatrix([]float64{2, 2, 3})
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("Before scale", S, "\n", row)
	S.ScaleByVec(S, row)
	fmt.Println("Scaled by row", S)
	col := row.T()
	fmt.Println("Transpose", col)
	S.ScaleByCol(S, col)
	fmt.Println("Scaled by col", S)
	rows2, _ := NewMatrix([]float64{2, 2, 2, 3, 3, 3})
	fmt.Println("Before adding 4", rows2)
	rows2.AddFloat(rows2, 4)
	fmt.Println("After adding 4", rows2)
}

func TestRowMod(Te *testing.T) {
	a := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}
	A, err := NewMatrix(a)
	if err != nil {
		Te.Error(err)
	}
	B := Zeros(5)
	B.DelVec(A, 3)
	fmt.Println("with and wihout row 3\n", A, "\n", B)
	fmt.Println("test for Unit")
	row, err := NewMatrix([]float64{2, 2, 3})
	if err != nil {
		Te.Error(err)
	}
	fmt.Println("Original vector", row)
	row.Unit(row)
	fmt.Println("Unitarized", row)

}

func TestEigen(Te *testing.T) {
	a := []float64{1, 2, 0, 2, 1, 0, 0, 0, 1}
	A, err := NewMatrix(a)
	if err != nil {
		Te.Error(err)
	}
	evecs, evals, err := EigenWrap(A, -1)
	fmt.Println(evecs, "\n", evals, "\n", err)
}
