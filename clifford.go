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

import "math"
import "runtime"

//import "fmt"

type paravector struct {
	Real  float64
	Imag  float64
	Vreal *VecMatrix
	Vimag *VecMatrix
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

//Takes a vector and creates a paravector. Uses copy so the vector is not affected
//By future changes to the paravector.
func paravectorFromVector(A,B *VecMatrix) *paravector {
	R := new(paravector)
	R.Real = 0 //I shouldnt need this
	R.Imag = 0
	R.Vreal = A
	R.Vimag = B
	return R
}

//Puts a copy fo the paravector P in the receiver
func (R *paravector) Copy(P *paravector) {
	R.Real = P.Real
	R.Imag = P.Imag
	R.Vreal.Copy(P.Vreal)
	R.Vimag.Copy(P.Vimag)
}

//Puts the reverse of paravector P in the received
func (R *paravector) reverse(P *paravector) {
	if P!=R{
		R.Copy(P)
	}
	R.Vimag.Scale(-1, R.Vimag)
}

//Puts  the normalized version of P in the receiver. If R and P are the same,
func (R *paravector) unit (P *paravector) {
	norm := 0.0
	norm += math.Pow(P.Real, 2) + math.Pow(P.Imag, 2)
	for i := 0; i < 3; i++ {
		norm += math.Pow(P.Vreal.At(0, i), 2) + math.Pow(P.Vimag.At(0, i), 2)
	}
	//fmt.Println("norm", norm)
	R.Real = P.Real / math.Sqrt(norm)
	R.Imag = P.Imag / math.Sqrt(norm)
	for i := 0; i < 3; i++ {
		R.Vreal.Set(0, i, P.Vreal.At(0, i)/math.Sqrt(norm))
		R.Vimag.Set(0, i, P.Vimag.At(0, i)/math.Sqrt(norm))
	}
	//fmt.Println("normalized", R)
}

//Clifford product of 2 paravectors, the imaginary parts are simply set to zero, since this is the case
//when rotating 3D real vectors. The proper Cliffor product is in fullCliProduct
func (R *paravector) cliProduct(A, B *paravector)  {
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
func (R *paravector) cliRotation(A, axis, tmp, tmp2 *paravector, angle float64) {
//	R := makeParavector()
	R.Real = math.Cos(angle / 2.0)
	for i := 0; i < 3; i++ {
		R.Vimag.Set(0, i, math.Sin(angle/2.0)*axis.Vreal.At(0, i))
	}
	R.reverse(R)
//	tmp:=makeParavector()
//	tmp2:=makeParavector()
	tmp.cliProduct(R, A)
	tmp2.cliProduct(tmp, R)
	R.Copy(tmp2)
}



//RotateSer takes the matrix Target and uses Clifford algebra to rotate each of its rows
//by angle radians around axis. Axis must be a 3D row vector. Target must be an N,3 matrix.
//The Ser in the name is from "serial". ToRot will be overwritten and returned
func RotateSer(Target,ToRot, axis *VecMatrix, angle float64) *VecMatrix {
	cake:=ZeroVecs(10) //Better ask for one chunk of memory than allocate 10 different times.
	R:=cake.VecView(0)
	Rrev:=cake.VecView(1)
	tmp:=cake.VecView(2)
	Rotated:=cake.VecView(3)
	itmp1:=cake.VecView(4)
	itmp2:=cake.VecView(5)
	itmp3:=cake.VecView(6)
	itmp4:=cake.VecView(7)
	itmp5:=cake.VecView(8)
	itmp6:=cake.VecView(9)
	RotateSerP(Target,ToRot,axis,tmp,R,Rrev,Rotated,itmp1,itmp2,itmp3,itmp4,itmp5,itmp6,angle)
	return ToRot
}


//RotateSerP takes the matrix Target and uses Clifford algebra to rotate each of its rows
//by angle radians around axis. Axis must be a 3D row vector. Target must be an N,3 matrix.
//The Ser in the name is from "serial". ToRot will be overwritten and returned. RotateSerP only allocates some floats but not
//any VecMatrix. Instead, it takes the needed intermediates as arguments, hence the "P" for "performance" If performance is not an issue,
//use RotateSer instead, it will perform the allocations for you and call this function. Notice that if you use this function directly
//you may have to zero at least some of the intermediates before reusing them.
func RotateSerP(Target,ToRot, axis,tmpv,Rv,Rvrev,Rotatedv, itmp1,itmp2,itmp3,itmp4,itmp5,itmp6 *VecMatrix, angle float64) {
	tarr, _ := Target.Dims()
	torotr := ToRot.NVecs()
	if tarr != torotr || Target.Dense == ToRot.Dense {
		panic("RotateSerP: Target and Res must have the same dimensions. Target and Res cannot reference the same matrix")
	}
	//Make the paravectors from the passed vectors:
	tmp:=paravectorFromVector(tmpv,itmp3)
	R:=paravectorFromVector(Rv,itmp4)
	Rrev:=paravectorFromVector(Rvrev,itmp5)
	Rotated:=paravectorFromVector(Rotatedv,itmp6)
	//That is with the building of temporary paravectors.
	paxis := paravectorFromVector(axis,itmp1)
	paxis.unit(paxis)
	R.Real = math.Cos(angle / 2.0)
	for i := 0; i < 3; i++ {
		R.Vimag.Set(0, i, math.Sin(angle/2.0)*paxis.Vreal.At(0, i))
	}
	Rrev.reverse(R)
	for i := 0; i < tarr; i++ {
		rowvec := Target.VecView(i)
		tmp.cliProduct(Rrev, paravectorFromVector(rowvec,itmp2))
		Rotated.cliProduct(tmp, R)
		ToRot.SetMatrix(i, 0, Rotated.Vreal)
	}
}


//Rotate takes the matrix Target and uses Clifford algebra to _concurrently_ rotate each
//of its rows by angle radians around axis. Axis must be a 3D row vector.
//Target must be an N,3 matrix. The result is returned.
func Rotate(Target, axis *VecMatrix, angle float64) *VecMatrix {
	Res := ZeroVecs(Target.NVecs())
	RotateP(Target, Res, axis, angle)
	return Res
}


//Rotate takes the matrix Target and uses Clifford algebra to _concurrently_ rotate each
//of its rows by angle radians around axis. Axis must be a 3D row vector.
//Target must be an N,3 matrix.

func RotateP(Target, Res, axis *VecMatrix, angle float64) {
	gorut := runtime.GOMAXPROCS(-1) //Do not change anything, only query
	rows := Target.NVecs()
	rrows := Res.NVecs()
	if rrows != rows || Target.Dense == Res.Dense {
		panic("RotateP: Target and Res must have the same dimensions. Target and Res cannot reference the same matrix.")
	}
	paxis := paravectorFromVector(axis,ZeroVecs(1))
	paxis.unit(paxis)
	R := makeParavector() //build the rotor (R)
	R.Real = math.Cos(angle / 2.0)
	for i := 0; i < 3; i++ {
		R.Vimag.Set(0, i, math.Sin(angle/2.0)*paxis.Vreal.At(0, i))
	}
	Rrev := makeParavector() //R-dagger
	Rrev.reverse(R)
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
			pv := paravectorFromVector(t2, t4)
			t3 := tmp3.VecView(i)
			for j := ini; j <= end; j++ {
				//Here we simply do R^dagger A R, and assign to the corresponding row.
				Rotated := paravectorFromVector(Res.VecView(j), t3)
				targetparavec := paravectorFromVector(Target.VecView(j), t1)
				pv.cliProduct(Rrev, targetparavec)
				Rotated.cliProduct(pv, R)
			}
			ended <- true
			return
		}(ini, end, i)
	}
	//Takes care of the concurrency
	for i := 0; i < gorut; i++ {
		<-ended
	}
	return 
}

/*
	rows, _ := Target.Dims()
	paxis := paravectorFromVector(axis,ZeroVecs(1))
	paxis.unit(paxis)
	R := makeParavector() //build the rotor (R)
	R.Real = math.Cos(angle / 2.0)
	for i := 0; i < 3; i++ {
		R.Vimag.Set(0, i, math.Sin(angle/2.0)*paxis.Vreal.At(0, i))
	}

	Rrev := R.reverse() // R-dagger
	Res := ZeroVecs(rows)
	ended := make(chan bool, rows)
	for i := 0; i < rows; i++ {
		go func(i int) {
			//Here we simply do R^dagger A R, and assign to the corresponding row.
			targetrow := Target.VecView(i)
			tmp := cliProduct(Rrev, paravectorFromVector(targetrow,ZeroVecs(1)))
			Rotated := cliProduct(tmp, R)
			//a,b:=Res.Dims() //debug
			//c,d:=Rotated.Vreal.Dims()
			//fmt.Println("rows",a,c,"cols",b,d,"i","rowss",)
			Res.SetMatrix(i, 0, Rotated.Vreal)
			ended <- true
			return
		}(i)
	}
	//Takes care of the concurrency
	for i := 0; i < rows; i++ {
		<-ended
	}
	return Res
}
*/

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
