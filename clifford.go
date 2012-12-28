/*
 * Clifford.go, part of gochem.
 * 
 * 
 * Copyright 2012 Janne Pesonen <janne.pesonen{at}helsinkiDOTfi> 
 * and Raul Mera <rmera{at}chemDOThelsinkiDOTfi> 
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
/***RM: Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/


	
package chem

import "fmt"
import  "github.com/skelterjohn/go.matrix"
import "math"

type paravector struct {
	Real float64
	Imag float64
	Vreal *matrix.DenseMatrix
	Vimag *matrix.DenseMatrix
	}

//creates a new paravector
func makeParavector() *paravector{
	R:=new(paravector)
	R.Real=0 //I shouldnt need this
	R.Imag=0
	R.Vreal=matrix.Zeros(1,3)
	R.Vimag=matrix.Zeros(1,3)
	return R
	}
//Takes a vector and creates a paravector. Uses copy so the vector is not affected
//By future changes to the paravector.
func paravectorFromVector(A *matrix.DenseMatrix) *paravector{
	R:=new(paravector)
	R.Real=0 //I shouldnt need this
	R.Imag=0
	R.Vreal=A.Copy()
	R.Vimag=matrix.Zeros(1,3)
	return R
}

//Returns a copy of the paravector.
func (P *paravector)Copy() *paravector{
	R:=new(paravector)
	R.Real=P.Real
	R.Imag=P.Imag
	R.Vreal=P.Vreal.Copy()
	R.Vimag=P.Vimag.Copy()
	return R
} 

//Returns the reverse of the paravector.
func (P *paravector)Reverse() *paravector{
	R:=new(paravector)
	R.Real=P.Real
	R.Imag=-1*P.Imag
	R.Vreal=P.Vreal.Copy()
	R.Vimag=P.Vimag.Copy()
	R.Vimag.Scale(-1)
	return R
} 

//Returns the normalized version of P.
func (P *paravector)Normalize() *paravector{
	R:=makeParavector()
	norm:=0.0
	norm+=math.Pow(P.Real,2) + math.Pow(P.Imag,2)
	for i:=0;i<3;i++{
		norm+=math.Pow(P.Vreal.Get(0,i),2)+math.Pow(P.Vimag.Get(0,i),2)
		}
	//fmt.Println("norm", norm)
	R.Real=P.Real/math.Sqrt(norm)
	R.Imag=P.Imag/math.Sqrt(norm)
	for i:=0;i<3;i++{
		R.Vreal.Set(0,i,P.Vreal.Get(0,i)/math.Sqrt(norm))
		R.Vimag.Set(0,i,P.Vimag.Get(0,i)/math.Sqrt(norm))
		}
	//fmt.Println("normalized", R)
	return R
}

//Clifford product of 2 paravectors, the imaginary parts are simply set to zero, since this is the case 
//when rotating 3D real vectors. The proper Cliffor product is in fullCliProduct
func cliProduct(A,B *paravector) *paravector{
	R:=makeParavector()
	R.Real=A.Real*B.Real-A.Imag*B.Imag
	for i:=0;i<3;i++{
		R.Real+=(A.Vreal.Get(0,i)*B.Vreal.Get(0,i) - A.Vimag.Get(0,i)*B.Vimag.Get(0,i)) 
	}
	R.Imag=A.Real*B.Imag + A.Imag*B.Real
	for i:=0;i<3;i++{
		R.Imag+=(A.Vreal.Get(0,i)*B.Vimag.Get(0,i) + A.Vimag.Get(0,i)*B.Vreal.Get(0,i)) 
	}
	//Now the vector part
	//First real
	R.Vreal.Set(0,0,A.Real*B.Vreal.Get(0,0) + B.Real*A.Vreal.Get(0,0) - A.Imag*B.Vimag.Get(0,0) - B.Imag*A.Vimag.Get(0,0) + 
	A.Vimag.Get(0,2)*B.Vreal.Get(0,1) - A.Vimag.Get(0,1)*B.Vreal.Get(0,2) + A.Vreal.Get(0,2)*B.Vimag.Get(0,1)- 
	A.Vreal.Get(0,1)*B.Vimag.Get(0,2))
	//Second real
	R.Vreal.Set(0,1,A.Real*B.Vreal.Get(0,1) + B.Real*A.Vreal.Get(0,1) - A.Imag*B.Vimag.Get(0,1) - B.Imag*A.Vimag.Get(0,1) +
	A.Vimag.Get(0,0)*B.Vreal.Get(0,2) - A.Vimag.Get(0,2)*B.Vreal.Get(0,0) + A.Vreal.Get(0,0)*B.Vimag.Get(0,2) - 
	A.Vreal.Get(0,2)*B.Vimag.Get(0,0))
	//Third real
	R.Vreal.Set(0,2,A.Real*B.Vreal.Get(0,2) + B.Real*A.Vreal.Get(0,2) - A.Imag*B.Vimag.Get(0,2) - B.Imag*A.Vimag.Get(0,2) +
	A.Vimag.Get(0,1)*B.Vreal.Get(0,0) - A.Vimag.Get(0,0)*B.Vreal.Get(0,1) + A.Vreal.Get(0,1)*B.Vimag.Get(0,0) - 
	A.Vreal.Get(0,0)*B.Vimag.Get(0,1))
	/*
	//First imag
	R.Vimag.Set(0,0,A.Real*B.Vimag.Get(0,0) + B.Real*A.Vimag.Get(0,0) + A.Imag*B.Vreal.Get(0,0) - B.Imag*A.Vreal.Get(0,0) +
	A.Vreal.Get(0,1)*B.Vreal.Get(0,2) - A.Vreal.Get(0,2)*B.Vreal.Get(0,1) + A.Vimag.Get(0,2)*B.Vimag.Get(0,1) -
	A.Vimag.Get(0,1)*B.Vimag.Get(0,2))
	//Second imag
	R.Vimag.Set(0,1,A.Real*B.Vimag.Get(0,1) + B.Real*A.Vimag.Get(0,1) + A.Imag*B.Vreal.Get(0,1) - B.Imag*A.Vreal.Get(0,1) +
	A.Vreal.Get(0,2)*B.Vreal.Get(0,0) - A.Vreal.Get(0,0)*B.Vreal.Get(0,2) + A.Vimag.Get(0,0)*B.Vimag.Get(0,2) - 
	A.Vimag.Get(0,2)*B.Vimag.Get(0,0))
	//Third imag
	R.Vimag.Set(0,2,A.Real*B.Vimag.Get(0,2) + B.Real*A.Vimag.Get(0,2) + A.Imag*B.Vreal.Get(0,2) - B.Imag*A.Vreal.Get(0,2) +
	A.Vreal.Get(0,0)*B.Vreal.Get(0,1) - A.Vreal.Get(0,1)*B.Vreal.Get(0,0) + A.Vimag.Get(0,1)*B.Vimag.Get(0,0) - 
	A.Vimag.Get(0,0)*B.Vimag.Get(0,1)) 
	*/
	//fmt.Println("R slido del horno", R)
	// A.Real, B.Vimag.Get(0,0), "g2", B.Real,A.Vimag.Get(0,0),"g3", A.Imag, B.Vreal.Get(0,0),"g4" ,B.Imag,A.Vreal.Get(0,0),
	//"g5", A.Vreal.Get(0,2), B.Vreal.Get(0,1), -1*A.Vreal.Get(0,1)*B.Vreal.Get(0,2), A.Vimag.Get(0,2)*B.Vimag.Get(0,1), -1*
    //A.Vimag.Get(0,1)*B.Vimag.Get(0,2))
	return R
	}
	
//cliRotation uses Clifford algebra to rotate a paravector Aby angle radians around axis. Returns the rotated 
//paravector. axis must be normalized.	
func cliRotation(A, axis *paravector, angle float64) *paravector{
	R:=makeParavector()
//	fmt.Println("Norm axis", axis) ////////////7
	R.Real=math.Cos(angle/2.0)
	for i:=0;i<3;i++{
		R.Vimag.Set(0,i,math.Sin(angle/2.0)*axis.Vreal.Get(0,i))
	}
	//fmt.Println("R", R)  ///////////////////////
	tmp:=cliProduct(R.Reverse(),A)
	Rotated:=cliProduct(tmp,R)
	//fmt.Println("rotated",Rotated) ////////////
	return Rotated
}

//CliRotate takes the matrix Target and uses Clifford algebra to rotate each of its rows 
//by angle radians around axis. Axis must be a 3D row vector. Target must be an N,3 matrix.
func CliRotate(Target, axis *matrix.DenseMatrix, angle float64) *matrix.DenseMatrix{
	fmt.Println("ax", axis)
	paxis:=paravectorFromVector(axis)
	fmt.Println("paxis:", paxis)
	paxis=paxis.Normalize()
	R:=matrix.Zeros(Target.Rows(),3)
	for i:=0;i<Target.Rows();i++{
	//	fmt.Println(i)
		tmp:=cliRotation(paravectorFromVector(Target.GetRowVector(i)),paxis,angle)
		R.SetMatrix(i,0,tmp.Vreal)
		}
	return R
	}

//CliRotate takes the matrix Target and uses Clifford algebra to _concurrently_ rotate each
//of its rows by angle radians around axis. Axis must be a 3D row vector. 
//Target must be an N,3 matrix.
func CliRotateConc(Target, axis *matrix.DenseMatrix, angle float64) *matrix.DenseMatrix{
	rows:=Target.Rows()
	fmt.Println("ax", axis)
	paxis:=paravectorFromVector(axis)
	fmt.Println("paxis:", paxis)
	paxis=paxis.Normalize()
	R:=matrix.Zeros(rows,3)
	ended:=make(chan bool,rows)
	for i:=0;i<rows;i++{
		go func(i int){
		tmp:=cliRotation(paravectorFromVector(Target.GetRowVector(i)),paxis,angle)
	//	fmt.Println("I from goruntine!", i)
		R.SetMatrix(i,0,tmp.Vreal) 
	//	R.Set(i,0,tmp.Vreal.Get(0,0))
	//	R.Set(i,1,tmp.Vreal.Get(0,1))
	//	R.Set(i,2,tmp.Vreal.Get(0,2))
		ended<-true
		return
		}(i)
	}
	for i:=0;i<rows;i++{
		<-ended
		}
	return R
}















//Clifford product of 2 paravectors.
func fullCliProduct(A,B *paravector) *paravector{
	R:=makeParavector()
	R.Real=A.Real*B.Real-A.Imag*B.Imag
	for i:=0;i<3;i++{
		R.Real+=(A.Vreal.Get(0,i)*B.Vreal.Get(0,i) - A.Vimag.Get(0,i)*B.Vimag.Get(0,i)) 
	}
	R.Imag=A.Real*B.Imag + A.Imag*B.Real
	for i:=0;i<3;i++{
		R.Imag+=(A.Vreal.Get(0,i)*B.Vimag.Get(0,i) + A.Vimag.Get(0,i)*B.Vreal.Get(0,i)) 
	}
	//Now the vector part
	//First real
	R.Vreal.Set(0,0,A.Real*B.Vreal.Get(0,0) + B.Real*A.Vreal.Get(0,0) - A.Imag*B.Vimag.Get(0,0) - B.Imag*A.Vimag.Get(0,0) + 
	A.Vimag.Get(0,2)*B.Vreal.Get(0,1) - A.Vimag.Get(0,1)*B.Vreal.Get(0,2) + A.Vreal.Get(0,2)*B.Vimag.Get(0,1)- 
	A.Vreal.Get(0,1)*B.Vimag.Get(0,2))
	//Second real
	R.Vreal.Set(0,1,A.Real*B.Vreal.Get(0,1) + B.Real*A.Vreal.Get(0,1) - A.Imag*B.Vimag.Get(0,1) - B.Imag*A.Vimag.Get(0,1) +
	A.Vimag.Get(0,0)*B.Vreal.Get(0,2) - A.Vimag.Get(0,2)*B.Vreal.Get(0,0) + A.Vreal.Get(0,0)*B.Vimag.Get(0,2) - 
	A.Vreal.Get(0,2)*B.Vimag.Get(0,0))
	//Third real
	R.Vreal.Set(0,2,A.Real*B.Vreal.Get(0,2) + B.Real*A.Vreal.Get(0,2) - A.Imag*B.Vimag.Get(0,2)- B.Imag*A.Vimag.Get(0,2)+
	A.Vimag.Get(0,1)*B.Vreal.Get(0,0) - A.Vimag.Get(0,0)*B.Vreal.Get(0,1) +  A.Vreal.Get(0,1)*B.Vimag.Get(0,0) - 
	A.Vreal.Get(0,0)*B.Vimag.Get(0,1))
	//First imag
	R.Vimag.Set(0,0,A.Real*B.Vimag.Get(0,0) + B.Real*A.Vimag.Get(0,0) + A.Imag*B.Vreal.Get(0,0) - B.Imag*A.Vreal.Get(0,0) +
	A.Vreal.Get(0,1)*B.Vreal.Get(0,2) - A.Vreal.Get(0,2)*B.Vreal.Get(0,1) + A.Vimag.Get(0,2)*B.Vimag.Get(0,1) -
	A.Vimag.Get(0,1)*B.Vimag.Get(0,2))
	//Second imag
	R.Vimag.Set(0,1,A.Real*B.Vimag.Get(0,1) + B.Real*A.Vimag.Get(0,1) + A.Imag*B.Vreal.Get(0,1) - B.Imag*A.Vreal.Get(0,1) +
	A.Vreal.Get(0,2)*B.Vreal.Get(0,0) - A.Vreal.Get(0,0)*B.Vreal.Get(0,2) + A.Vimag.Get(0,0)*B.Vimag.Get(0,2) - 
	A.Vimag.Get(0,2)*B.Vimag.Get(0,0))
	//Third imag
	R.Vimag.Set(0,2,A.Real*B.Vimag.Get(0,2) + B.Real*A.Vimag.Get(0,2) + A.Imag*B.Vreal.Get(0,2) - B.Imag*A.Vreal.Get(0,2) +
	A.Vreal.Get(0,0)*B.Vreal.Get(0,1) - A.Vreal.Get(0,1)*B.Vreal.Get(0,0) + A.Vimag.Get(0,1)*B.Vimag.Get(0,0) - 
	A.Vimag.Get(0,0)*B.Vimag.Get(0,1))
	//fmt.Println("R slido del horno", R)
	// A.Real, B.Vimag.Get(0,0), "g2", B.Real,A.Vimag.Get(0,0),"g3", A.Imag, B.Vreal.Get(0,0),"g4" ,B.Imag,A.Vreal.Get(0,0),
	//"g5", A.Vreal.Get(0,2), B.Vreal.Get(0,1), -1*A.Vreal.Get(0,1)*B.Vreal.Get(0,2), A.Vimag.Get(0,2)*B.Vimag.Get(0,1), -1*
    //A.Vimag.Get(0,1)*B.Vimag.Get(0,2))
    
	return R
	}
