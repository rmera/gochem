/*
 * untitled.go
 * 
 * Copyright 2012 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
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
 */


	
package chem

import "fmt"
import  "github.com/skelterjohn/go.matrix"

type paravector struct {
	Real float64
	Imag float64
	Vreal *matrix.DenseMatrix
	Vimag *matrix.DenseMatrix
	}




func CliProduct(A,B *paravector) paravector{
	var  R paravector
	R.Real:=A.Real*B.Real-A.Imag*B.Imag
	for i:=0;i<3;i++{
		R.Real+=(A.Vreal.Get(0,i)*B.Vreal.Get(0,i) - A.Vimag.Get(0,i)*B.Vimag.Get(0,i)) 
		}
	R.Imag=A.Real*B.Imag + A.Imag*B.Real
	for i:=0;i<3;i++{
		R.Imag+=(A.Vreal.Get(0,i)*B.Vimag.Get(0,i) + A.Vimag.Get(0,i)*B.Vreal.Get(0,i)) 
		}
	//Now the vector part
	//first real
	real0:=A.Real*B.Vreal.Get(0,0) + B.Real*A.Vreal.Get(0,0) - A.Imag*B.Vimag.Get(0,0) - A.Imag*B.Vimag.Get(0,0) 
	+ A.Vimag.Get(0,2)*B.Vreal(0,1) - A.Vimag.Get(0,1)*B.Vreal(0,2) + A.Vreal.Get(0,2)*B.Vimag(0,1)
	- A.Vreal.Get(0,1)*B.Vimag(0,2)
	R.Vreal.Set(0,0,real0)
	//second real
	}
	
	
	
