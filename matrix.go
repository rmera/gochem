/*
 * matrix.go, part of gochem.
 * 
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
 * 
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 * 
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

	
package chem



import "github.com/skelterjohn/go.matrix"
import "math"
import "fmt"

//These ones are basic math, belonging more to the go.matrix package
//If there is something similar already made
//in go.matrix this functions will be deleted. Otherwise they could be
//made methods for DenseMatrix and included in go.matrix

//Some of this functions don't return error messages because they are meant to
//Be inserted in mathematical expressions and thus they need to return only one value.

/*SetRows replaces the coordinates of atoms in the positions given by atomlist with the coordinates
in newcoords (in order) If atomlist contains a single element, it replaces as many coordinates
as given in newcoords, starting 
at the element in atomlist. In the latter case, the function checks that there are enough coordinates to
replace and panics if not. */
func SetRows(target, newcoords *matrix.DenseMatrix, atomlist []int) *matrix.DenseMatrix{
	//If supplies a list with one number, the newcoords will replace the old coords
	//Starting that number. We do check that you don't put more coords than spaces we have.
	returned:=target.Copy()
	if len(atomlist)==1{
		if newcoords.Rows()>target.Rows()-atomlist[0]-1{
			panic(fmt.Sprintf("Cant replace starting from position %d: Not enough atoms in molecule", atomlist[0]))
			} 
		returned.SetMatrix(atomlist[0],0,newcoords)
		return returned
		}
	//If the list has more than one atom
	lenatoms:=target.Rows()	
	for k,j:=range(atomlist){
		if j>lenatoms-1{
			panic(fmt.Sprintf("Requested position number: %d (%d) out of range",k,j))
			}
		returned.SetMatrix(j,0,newcoords.GetRowVector(k))
		}
	return returned
	}
	



//Somerows, given a list of ints and the desired frame, returns an slice matrix.DenseMatrix
//containing the coordinates of the atoms with the corresponding index.
//This function returns a copy, not a reference, so changes to the returned matrix
//don't alter the original. It check for correctness of the frame and the
//Atoms requested.
func SomeRows(M *matrix.DenseMatrix, clist []int) (*matrix.DenseMatrix){
	rows:=M.Rows()
	ret:=make([][]float64,0,len(clist))
	for k,j:=range(clist){
		if j>rows-1{
			panic(fmt.Sprintf("Coordinate requested (Number: %d, value: %d) out of range!",k,j))
			}
		tmp:=M.GetRowVector(j).Array()
		if len(tmp)!=3{
			panic(fmt.Sprintf("Coordinate %d has %d components instead of 3",k, len(tmp)))
			}
		ret=append(ret,tmp)
		}
	return matrix.MakeDenseMatrixStacked(ret)
	}

//Unitarize takes a vector and divides it by its norm
//thus obtaining an unitary vector pointing in the same direction as 
//vector.
func Unitarize(vector *matrix.DenseMatrix) *matrix.DenseMatrix{
	norm:=vector.TwoNorm()
	norm=1.0/norm
	F:=vector.Copy()
	F.Scale(norm)
	return F
	}

//AddRow adds the row vector row to each row of the matrix big, in place. Both need the same ammount of columns.
func AddRow(big,row *matrix.DenseMatrix)(error){
	bigrows:=big.Rows()
	if big.Cols() != row.Cols() || row.Rows()!=1{
		return fmt.Errorf("Ill-formed matrices for multiplication")
		}
	for i:=0;i<bigrows;i++{
		j:=big.GetRowVector(i)
		j.Add(row)
		}
	return nil
	}

//SubRow substract the row bector row from each row of the matrix big, in place. Bothe ned the same ammount 
//of coulmn. This is not an efficient function since it contains two function calls, one to the DenseMatrix
//method Scale() and other to the function AddRow().
func SubRow(big,row *matrix.DenseMatrix)(error){
	row2:=row.Copy()
	row2.Scale(-1)
	err:=AddRow(big,row2)
	return err  //nil if no error, something else otherwise
	}

//Cross3D Takes 2 3-len column or row vectors and returns a column or a row
//vector, respectively, with the Cross product of them.
func Cross3D(a,b *matrix.DenseMatrix)(*matrix.DenseMatrix,error){
	ac:=a.Cols()
	ar:=a.Rows()
	bc:=b.Cols()
	br:=b.Rows()
	if ac != bc || ar != br {
		return nil, fmt.Errorf("Malformed vectors for cross product")
		}
	if ac!=3 {
		//Ok, Im sure one can do this better.
		c:=a.Transpose()
		d:=b.Transpose()
		e:=Cross3DRow(c,d)
		f:=e.Transpose()
		return f, nil
		}
	if ar!=3 {
		return nil, fmt.Errorf("Malformed vectors for cross product")
		}
	return Cross3DRow(a,b), nil
	}

//Cross3DRow returns the cross product of 2 row vectors. No error checking!
func Cross3DRow(a,b *matrix.DenseMatrix)*matrix.DenseMatrix{
	vec:=make([]float64,3,3)
	vec[0]=a.Get(0,1)*b.Get(0,2) - a.Get(0,2)*b.Get(0,1)
	vec[1]=a.Get(0,2)*b.Get(0,0) - a.Get(0,0)*b.Get(0,2)
	vec[2]=a.Get(0,0)*b.Get(0,1) - a.Get(0,1)*b.Get(0,0)
	return matrix.MakeDenseMatrix(vec,1,3)
	}


//InvSqrt return the inverse of the square root of val, or zero if
//|val|<appzero. Returns -1 if val is negative 
func InvSqrt(val float64) float64{
	if math.Abs(val)<=appzero{  //Not sure if need the 
		return 0
		}else if val<0{  //negative
		return  -1  //might change
		}
	return 1.0/math.Sqrt(val)	
	}

//KronekerDelta is a naive implementation of the kroneker delta function.
func KronekerDelta(a,b float64) float64{
	if math.Abs(a-b)<=appzero{
		return 1
		}
	return 0	
	}
	
	

//Dot returns the dot product between 2 vectors or matrices. Just the sum of the 
//Element-wise multiplication. In this case returning error
//makes it problematic to use in complex operations, so this returns -1
//when problems.
func Dot(A, B *matrix.DenseMatrix) float64{
	var err error
	if A.Cols()!=B.Cols() || A.Rows()!=B.Rows(){
		panic("Matrices must have the same dimmension to obtain dot product")
		}
	//For some crazy reason if the F variable is called C, I get a
	//"ScaleMatrixDense undeclared" error at compile time :S
	F := A.Copy()
	err=F.ScaleMatrixDense(B)
	if err!=nil{
		panic(err.Error())
		}
	return DMSummation(F)
	}

//DMPowInPlace raises the DenseMatrix A, element-wise, to the nth power.  
func DMPowInPlace(A *matrix.DenseMatrix, n float64){
	for i:=0;i<A.Rows();i++{
		for j:=0;j<A.Cols();j++{
		//	fmt.Println(i,j, A.Rows(),A.Cols())
			A.Set(i,j,math.Pow(A.Get(i,j),n))
			}
		}
	}

//DMPow returns B^n. It does not modify B.
func DMPow(B *matrix.DenseMatrix, n float64) *matrix.DenseMatrix{
	A:=B.Copy()
	for i:=0;i<A.Rows();i++{
		for j:=0;j<A.Cols();j++{
		//	fmt.Println(i,j, A.Rows(),A.Cols())
			A.Set(i,j,math.Pow(A.Get(i,j),n))
			}
		}
	return A
	}

//DMScaleByCol scales each column of matrix A by Col.
func DMScaleByCol(A, Col *matrix.DenseMatrix) error{
	Rows:=A.Rows()
	if Rows != Col.Rows() || Col.Cols()>1{
		return fmt.Errorf("Malformed matrices for scaling")
		}
	for i:=0;i<Rows;i++{
		A.ScaleRow(i,Col.Get(i,0))
		}
	return nil
	}

//DMScaleByRow each row of matrix A by Row.
func DMScaleByRow(A, Row *matrix.DenseMatrix) error{
	Cols:=A.Cols()
	if Cols != Row.Cols() || Row.Rows()>1{
		return fmt.Errorf("Malformed matrices for scaling")
		}
	for i:=0;i<Cols;i++{
		mult:=Row.Get(0,i)
		for j:=0;j<A.Rows();j++{
			A.Set(j,i,mult*A.Get(j,i))
			}
		}
	return nil
	}

//DMSummation returns the sum of all elements in matrix A.
func DMSummation(A *matrix.DenseMatrix) float64 {
	Rows:=A.Rows()
	Cols:=A.Cols()
	var sum float64
	for i:=0;i<Cols;i++{
		for j:=0;j<Rows;j++{
			sum+=A.Get(j,i)
			}
		}
	return sum
	}
	
//DMaddfloat returns a matrix which elements are those of matrix A plus the float B.
func DMaddfloat(A *matrix.DenseMatrix, B float64)*matrix.DenseMatrix{
	Rows:=A.Rows()
	Cols:=A.Cols()
	copy:=matrix.MakeDenseCopy(A)
	for i:=0;i<Cols;i++{
		for j:=0;j<Rows;j++{
			copy.Set(j,i,(A.Get(j,i)+B))
			}
		}
	return copy
	}


//DMDelRow removes Row i from matrix A. 
func DMDelRow(A *matrix.DenseMatrix,i int)(*matrix.DenseMatrix, error){
		//Im not sure if copying takes place.
		var err error
		lenght:=A.Rows()
		upper:=A.GetMatrix(0,0,i,3)
		lower:=A.GetMatrix(i+1,0,lenght-i-1,3)
		A,err=upper.Stack(lower)
		return A, err
	}







