// +build matrix

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

package chem

import "testing"
import "fmt"

func TestGeo(Te *testing.T) {
	a := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9}
	A := NewCoords(a)
	ar,_:= A.Dims()
	T := ZeroVecs(ar)
	T.T(A)
	B := gnEye(ar)
	//B.Clone(A)
	T.Mul(A, B)
	E := ZeroVecs(ar)
	E.MulElem(A, B)
	fmt.Println(T, "\n", T, "\n", A, "\n", B, "\n", ar, ac, A.Sum())
	View := ZeroVecs(1)
	View.VecView(A, 0)
	View.Set(0, 0, 100)
	fmt.Println("View\n", A, "\n", View)

}

func TestSomeVecs(Te *testing.T) {
	a := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}
	A := NewCoords(a, 6, 3)
	B := Zeros(3, 3) //We should cause an error by modifying this.
	cind := []int{1, 3, 5}
	err := B.SomeVecsSafe(A, cind)
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
	A := NewCoords(a)
	B := Zeros(6, 3)
	A.Scale(3, A)
	B.Scale(2, A)
	fmt.Println(A, "\n", B)
	Row := NewCoords([]float64{10, 20, 30})
	A.AddRow(A, Row)
	fmt.Println("Additions")
	fmt.Println(A)
	A.SubRow(A, Row)
	fmt.Println(A)
	B.Pow(A, 2)
	fmt.Println("Squared", A, "\n", B)
	b := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9}
	S := NewCoords(b)
	row := NewCoords([]float64{2, 2, 3})
	fmt.Println("Before scale", S, "\n", row)
	S.ScaleByRow(S, row)
	fmt.Println("Scaled by row", S)
	col := Zeros(3, 1)
	col.T(row)
	fmt.Println("Transpose", col)
	S.ScaleByCol(S, col)
	fmt.Println("Scaled by col", S)
	rows2 := NewCoords([]float64{2, 2, 2, 3, 3, 3})
	fmt.Println("Before adding 4", rows2)
	rows2.AddFloat(rows2, 4)
	fmt.Println("After adding 4", rows2)
}

func TestRowMod(Te *testing.T) {
	a := []float64{1.0, 2.0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}
	A := NewCoords(a)
	B := Zeros(5, 3)
	B.DelRow(A, 3)
	fmt.Println("with and wihout row 3\n", A, "\n", B)
	fmt.Println("test for Unit")
	row := NewCoords([]float64{2, 2, 3})
	fmt.Println("Original vector", row)
	row.Unit(row)
	fmt.Println("Unitarized", row)

}

func TestEigen(Te *testing.T) {
	a := []float64{1, 2, 0, 2, 1, 0, 0, 0, 1}
	A := NewCoords(a)
	evecs, evals, err := gnEigen(A, -1)
	fmt.Println(evecs, "\n", evals, "\n", err)
	U, s, V, err := gnSVD(A)
	fmt.Println("U\n", U, "\nsigma\n", s, "V\n", V, "\n", err)

}
