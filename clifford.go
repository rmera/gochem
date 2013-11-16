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

import (
	"math"
	"runtime"
	"fmt"
)

//import "fmt"

//type dvector struct{}
//func (d *dvector)At(a,b int) float64{
//	return 0.0
//}

type paravector struct {
	Real  float64
	Imag  float64
	Vreal *VecMatrix
	Vimag *VecMatrix //cheating
}

//creates a new paravector
func makeParavector() *paravector {
	R := new(paravector)
	R.Real = 0 //I shouldnt need this
	R.Imag = 0
	R.Vreal = ZeroVecs(1)
	R.Vimag = ZeroVecs(1)
	return R
}

//Takes a vector and creates a paravector.
func paravectorFromVector(A *VecMatrix) *paravector {
	R := new(paravector)
	R.Real = 0 //I shouldnt need this
	R.Imag = 0
	R.Vreal = A
	R.Vimag = ZeroVecs(1)
	return R
}

//Takes 2 vectors and creates a paravector. The second vector becomes the imaginary part.
func paravectorFromVectors(A, B *VecMatrix) *paravector {
	R := new(paravector)
	R.Real = 0 //I shouldnt need this
	R.Imag = 0
	R.Vreal = A
	R.Vimag = B
	return R
}

//Copies the given paravector into the receiver
func (R *paravector) Copy(P *paravector) {
	R.Real = P.Real
	R.Imag = P.Imag
	//	R.Vreal = ZeroVecs(1)
	//	R.Vimag = ZeroVecs(1)
	R.Vreal.Copy(P.Vreal)
	R.Vimag.Copy(P.Vimag)
	//	return R
}

//Puts a reversed paravector into the received
func (R *paravector) reverse(P *paravector) {
	R.Copy(P)
	R.Vimag.Scale(-1, R.Vimag)
	//	return R
}

//Puts a normalized version of P in the receiver.
func (R *paravector) unit(P *paravector) {
	norm := 0.0
	norm += math.Pow(P.Real, 2) + math.Pow(P.Imag, 2)
	for i := 0; i < 3; i++ {
		norm += math.Pow(P.Vreal.At(0, i), 2) + math.Pow(P.Vimag.At(0, i), 2)
	}
	//fmt.Println("norm", norm)
	N := 1 / math.Sqrt(norm)
	R.Real *= N
	R.Imag *= N
	R.Vreal.Scale(N, R.Vreal)
	R.Vimag.Scale(N, R.Vimag)
	/*
		for i := 0; i < 3; i++ {
			R.Vreal.Set(0, i, P.Vreal.At(0, i)/math.Sqrt(norm))
			R.Vimag.Set(0, i, P.Vimag.At(0, i)/math.Sqrt(norm))
		}
		//fmt.Println("normalized", R)
		return R
	*/
}

//Clifford product of 2 paravectors, the imaginary parts are simply set to zero, since this is the case
//when rotating 3D real vectors. The proper Cliffor product is in fullCliProduct
func (R *paravector) cliProduct(A, B *paravector) {
	R.Real = A.Real*B.Real - A.Imag*B.Imag
	Ar := []float64{A.Vreal.At(0, 0), A.Vreal.At(0, 1), A.Vreal.At(0, 2)}
	Ai := []float64{A.Vimag.At(0, 0), A.Vimag.At(0, 1), B.Vimag.At(0, 2)}

	Br := []float64{B.Vreal.At(0, 0), B.Vreal.At(0, 1), B.Vreal.At(0, 2)}
	Bi := []float64{B.Vimag.At(0, 0), B.Vimag.At(0, 1), B.Vimag.At(0, 2)}

	for i := 0; i < 3; i++ {
		R.Real += (Ar[i]*Br[i] - Ai[i]*Bi[i])
	}
	R.Imag = A.Real*B.Imag + A.Imag*B.Real
	for i := 0; i < 3; i++ {
		R.Imag += (Ar[i]*Bi[i] + Ai[i]*Br[i])
	}
	//Now the vector part
	//First real
	R.Vreal.Set(0, 0, A.Real*Br[0]+B.Real*Ar[0]-A.Imag*Bi[0]-B.Imag*Ai[0]+
		Ai[2]*Br[1]-Ai[1]*Br[2]+Ar[2]*Bi[1]-
		Ar[1]*Bi[2])
	//Second real
	R.Vreal.Set(0, 1, A.Real*Br[1]+B.Real*Ar[1]-A.Imag*Bi[1]-B.Imag*Ai[1]+
		Ai[0]*Br[2]-Ai[2]*Br[0]+Ar[0]*Bi[2]-
		Ar[2]*Bi[0])
	//Third real
	R.Vreal.Set(0, 2, A.Real*Br[2]+B.Real*Ar[2]-A.Imag*Bi[2]-B.Imag*Ai[2]+
		Ai[1]*Br[0]-Ai[0]*Br[1]+Ar[1]*Bi[0]-
		Ar[0]*Bi[1])
	/*
		//First imag
		R.Vimag.Set(0,0,A.Real*B.Vimag.At(0,0) + B.Real*A.Vimag.At(0,0) + A.Imag*B.Vreal.At(0,0) - B.Imag*A.Vreal.At(0,0) +
		A.Vreal.At(0,1)*B.Vreal.At(0,2) - A.Vreal.At(0,2)*B.Vreal.At(0,1) + A.Vimag.At(0,2)*B.Vimag.At(0,1) -
		A.Vimag.At(0,1)*B.Vimag.At(0,2))
		//Second imag
		R.Vimag.Set(0,1,A.Real*B.Vimag.At(0,1) + B.Real*A.Vimag.At(0,1) + A.Imag*B.Vreal.At(0,1) - B.Imag*A.Vreal.At(0,1) +
		A.Vreal.At(0,2)*B.Vreal.At(0,0) - A.Vreal.At(0,0)*B.Vreal.At(0,2) + A.Vimag.At(0,0)*B.Vimag.At(0,2) -
		A.Vimag.At(0,2)*B.Vimag.At(0,0))
		//Third imag
		R.Vimag.Set(0,2,A.Real*B.Vimag.At(0,2) + B.Real*A.Vimag.At(0,2) + A.Imag*B.Vreal.At(0,2) - B.Imag*A.Vreal.At(0,2) +
		A.Vreal.At(0,0)*B.Vreal.At(0,1) - A.Vreal.At(0,1)*B.Vreal.At(0,0) + A.Vimag.At(0,1)*B.Vimag.At(0,0) -
		A.Vimag.At(0,0)*B.Vimag.At(0,1))
	*/
	//fmt.Println("R slido del horno", R)
	// A.Real, B.Vimag.At(0,0), "g2", B.Real,A.Vimag.At(0,0),"g3", A.Imag, B.Vreal.At(0,0),"g4" ,B.Imag,A.Vreal.At(0,0),
	//"g5", A.Vreal.At(0,2), B.Vreal.At(0,1), -1*A.Vreal.At(0,1)*B.Vreal.At(0,2), A.Vimag.At(0,2)*B.Vimag.At(0,1), -1*
	//A.Vimag.At(0,1)*B.Vimag.At(0,2))
	//	return R
}

//cliRotation uses Clifford algebra to rotate a paravector Aby angle radians around axis. Returns the rotated
//paravector. axis must be normalized.
func (R *paravector) cliRotation(A, axis, tmp1, tmp2 *paravector, angle float64) {
	angle /= 2.0
	R.Real = math.Cos(angle)
	for i := 0; i < 3; i++ {
		R.Vimag.Set(0, i, math.Sin(angle)*axis.Vreal.At(0, i))
	}
	tmp1.reverse(R)
	tmp2.cliProduct(tmp1, A)
	R.cliProduct(tmp2, R)
}

/*
//RotateSer takes the matrix Target and uses Clifford algebra to rotate each of its rows
//by angle radians around axis. Axis must be a 3D row vector. Target must be an N,3 matrix.
//The Ser in the name is from "serial". The CliRotate is concurrent.
func RotateSer(Target, axis *VecMatrix, angle float64) *VecMatrix {
	tarr, _ := Target.Dims()
	paxis := paravectorFromVector(axis)
	paxis = paxis.unit()
	R := makeParavector()
	R.Real = math.Cos(angle / 2.0)
	for i := 0; i < 3; i++ {
		R.Vimag.Set(0, i, math.Sin(angle/2.0)*paxis.Vreal.At(0, i))
	}
	Rrev := R.reverse()
	Res := ZeroVecs(tarr)
	for i := 0; i < tarr; i++ {
		rowvec := Target.VecView(i)
		Rotated:=paravectorFromVectors(Res.VecView(i),t1)
		t2.cliProduct(Rrev, paravectorFromVector(rowvec,t2))
		Rotated := cliProduct(t2, R)
	}
	return Res
}
*/

//Rotate takes the matrix Target and uses Clifford algebra to _concurrently_ rotate each
//of its rows by angle radians around axis. Axis must be a 3D row vector.
//Target must be an N,3 matrix.
func Rotate(Target, axis *VecMatrix, angle float64) *VecMatrix {
	gorut := runtime.GOMAXPROCS(-1) //Do not change anything, only query
	rows, _ := Target.Dims()
	paxis := paravectorFromVector(axis)
	paxis.unit(paxis)
	R := makeParavector() //build the rotor (R)
	R.Real = math.Cos(angle / 2.0)
	for i := 0; i < 3; i++ {
		R.Vimag.Set(0, i, math.Sin(angle/2.0)*paxis.Vreal.At(0, i))
	}
	Rrev := makeParavector() //R-dagger
	Rrev.reverse(R)
	Res := ZeroVecs(rows) //The returned matrix
	ended := make(chan bool, gorut)
	//Just temporal skit to be used by the gorutines
	tmp1 := ZeroVecs(gorut)
	tmp2 := ZeroVecs(gorut)
	tmp3 := ZeroVecs(gorut)
	tmp4 := ZeroVecs(gorut)
	fragmentlen := int(math.Ceil(float64(rows) / float64(gorut))) //len of the fragment of target that each gorutine will handle
	for i := 0; i < gorut; i++ {
		//These are the limits of the fragment of Target in which the gorutine will operate
		ini := i * fragmentlen
		end := i*fragmentlen + (fragmentlen - 1)
		if i == gorut-1 {
			end = rows - 1 //The last fragment may be smaller than fragmentlen
		}
		go func(ini, end, i int) {
			t1 := tmp1.VecView(i)
			t2 := tmp2.VecView(i)
			t4 := tmp4.VecView(i)
			pv := paravectorFromVectors(t2, t4)
			t3 := tmp3.VecView(i)
			for j := ini; j <= end; j++ {
				//Here we simply do R^dagger A R, and assign to the corresponding row.
				Rotated := paravectorFromVectors(Res.VecView(j), t3)
				targetparavec := paravectorFromVectors(Target.VecView(j),t1)
				pv.cliProduct(Rrev, targetparavec)
				Rotated.cliProduct(pv, R)
				//a,b:=Res.Dims() //debug
				//c,d:=Rotated.Vreal.Dims()
				//fmt.Println("rows",a,c,"cols",b,d,"i","rowss",)
		//		t1.Set(0, 0, 0.0) //This could be more expensive than simpy getting new vectors. Profile!
		//		t1.Set(0, 1, 0.0)
		//		t1.Set(0, 2, 0.0)
			}
			ended <- true
			return
		}(ini, end, i)
	}
	//Takes care of the concurrency
	for i := 0; i < gorut; i++ {
		<-ended
	}
	fmt.Println("gorut",gorut)
	return Res
}

//Clifford product of 2 paravectors.
func fullCliProduct(A, B *paravector) *paravector {
	R := makeParavector()
	R.Real = A.Real*B.Real - A.Imag*B.Imag
	for i := 0; i < 3; i++ {
		R.Real += (A.Vreal.At(0, i)*B.Vreal.At(0, i) - A.Vimag.At(0, i)*B.Vimag.At(0, i))
	}
	R.Imag = A.Real*B.Imag + A.Imag*B.Real
	for i := 0; i < 3; i++ {
		R.Imag += (A.Vreal.At(0, i)*B.Vimag.At(0, i) + A.Vimag.At(0, i)*B.Vreal.At(0, i))
	}
	//Now the vector part
	//First real
	R.Vreal.Set(0, 0, A.Real*B.Vreal.At(0, 0)+B.Real*A.Vreal.At(0, 0)-A.Imag*B.Vimag.At(0, 0)-B.Imag*A.Vimag.At(0, 0)+
		A.Vimag.At(0, 2)*B.Vreal.At(0, 1)-A.Vimag.At(0, 1)*B.Vreal.At(0, 2)+A.Vreal.At(0, 2)*B.Vimag.At(0, 1)-
		A.Vreal.At(0, 1)*B.Vimag.At(0, 2))
	//Second real
	R.Vreal.Set(0, 1, A.Real*B.Vreal.At(0, 1)+B.Real*A.Vreal.At(0, 1)-A.Imag*B.Vimag.At(0, 1)-B.Imag*A.Vimag.At(0, 1)+
		A.Vimag.At(0, 0)*B.Vreal.At(0, 2)-A.Vimag.At(0, 2)*B.Vreal.At(0, 0)+A.Vreal.At(0, 0)*B.Vimag.At(0, 2)-
		A.Vreal.At(0, 2)*B.Vimag.At(0, 0))
	//Third real
	R.Vreal.Set(0, 2, A.Real*B.Vreal.At(0, 2)+B.Real*A.Vreal.At(0, 2)-A.Imag*B.Vimag.At(0, 2)-B.Imag*A.Vimag.At(0, 2)+
		A.Vimag.At(0, 1)*B.Vreal.At(0, 0)-A.Vimag.At(0, 0)*B.Vreal.At(0, 1)+A.Vreal.At(0, 1)*B.Vimag.At(0, 0)-
		A.Vreal.At(0, 0)*B.Vimag.At(0, 1))
	//First imag
	R.Vimag.Set(0, 0, A.Real*B.Vimag.At(0, 0)+B.Real*A.Vimag.At(0, 0)+A.Imag*B.Vreal.At(0, 0)-B.Imag*A.Vreal.At(0, 0)+
		A.Vreal.At(0, 1)*B.Vreal.At(0, 2)-A.Vreal.At(0, 2)*B.Vreal.At(0, 1)+A.Vimag.At(0, 2)*B.Vimag.At(0, 1)-
		A.Vimag.At(0, 1)*B.Vimag.At(0, 2))
	//Second imag
	R.Vimag.Set(0, 1, A.Real*B.Vimag.At(0, 1)+B.Real*A.Vimag.At(0, 1)+A.Imag*B.Vreal.At(0, 1)-B.Imag*A.Vreal.At(0, 1)+
		A.Vreal.At(0, 2)*B.Vreal.At(0, 0)-A.Vreal.At(0, 0)*B.Vreal.At(0, 2)+A.Vimag.At(0, 0)*B.Vimag.At(0, 2)-
		A.Vimag.At(0, 2)*B.Vimag.At(0, 0))
	//Third imag
	R.Vimag.Set(0, 2, A.Real*B.Vimag.At(0, 2)+B.Real*A.Vimag.At(0, 2)+A.Imag*B.Vreal.At(0, 2)-B.Imag*A.Vreal.At(0, 2)+
		A.Vreal.At(0, 0)*B.Vreal.At(0, 1)-A.Vreal.At(0, 1)*B.Vreal.At(0, 0)+A.Vimag.At(0, 1)*B.Vimag.At(0, 0)-
		A.Vimag.At(0, 0)*B.Vimag.At(0, 1))
	//fmt.Println("R slido del horno", R)
	// A.Real, B.Vimag.At(0,0), "g2", B.Real,A.Vimag.At(0,0),"g3", A.Imag, B.Vreal.At(0,0),"g4" ,B.Imag,A.Vreal.At(0,0),
	//"g5", A.Vreal.At(0,2), B.Vreal.At(0,1), -1*A.Vreal.At(0,1)*B.Vreal.At(0,2), A.Vimag.At(0,2)*B.Vimag.At(0,1), -1*
	//A.Vimag.At(0,1)*B.Vimag.At(0,2))

	return R
}
