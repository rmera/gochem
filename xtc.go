// +build gromacs

/*
 * untitled.go
 * 
 * Copyright 2012 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License  as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */
  /* 
 * 
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 * 
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

	
package chem

// #cgo CFLAGS:  -I.
// #cgo LDFLAGS:  -L.  -lnsl -lm   -lxdrfile
//#include <xtc.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <xdrfile_xtc.h>
//#include <xdrfile_trr.h>
//#include <xdrfile.h>
import "C"
import "fmt"
import "github.com/skelterjohn/go.matrix"

func ReadXtcFrames(ini, end, skip int, filename string)([]*matrix.DenseMatrix,int, error) {
	Coords:=make([]*matrix.DenseMatrix,0,1) // I might attempt to give the capacity later
	natoms,err:=XtcCountAtoms(filename)
	if err!=nil{
		return nil, 0, err
		}
	fp,_:=XtcOpen(filename) //We already tested that the name is ok, no need to catch this error
	ccoords:=make([]C.float,natoms*3)
	i:=0
	for ;;i++{
		if i>end {
			break
			}
		if i<ini || (i-(1+ini))%skip!=0{
			err:=xtcGetFrameEfficientDrop(fp,ccoords,natoms)
			if err!=nil{
				if err.Error()=="No more frames"{
					break //No more frames is not really an error
					}else{
					return Coords, i, err	
					}
				}
			}else{
			coords,err:=xtcGetFrameEfficient(fp,ccoords,natoms)	
			if err!=nil{
				if err.Error()=="No more frames"{
					break //No more frames is not really an error
					}else{
					return Coords, i, err
					}
				}
			Coords=append(Coords,matrix.MakeDenseMatrix(coords,natoms,3))
			}
		}
	return Coords, i, nil
	}
	



func XtcCountAtoms(name string)(int, error){
	Cfilename:=C.CString(name)
	Cnatoms:=C.read_natoms(Cfilename)
	natoms:=int(Cnatoms)
	return natoms, nil
	}
	
func  XtcOpen(name string)(*C.XDRFILE, error){
	Cfilename:=C.CString(name)
	fp:=C.openfile(Cfilename)
	if fp==nil{
		return nil, fmt.Errorf("Unable to open xtc file")
		}
	return fp, nil
	}

func XtcGetFrame(fp *C.XDRFILE,natoms int)([]float64,error){
	totalcoords:=natoms*3
	cnatoms:=C.int(natoms)
	Ccoords:= make([]C.float,totalcoords)
	worked:=C.get_coords(fp,&Ccoords[0],cnatoms)
	if worked==11{
		goCoords:=make([]float64,1)
		return goCoords, fmt.Errorf("No more frames")
		}
	if worked!=0{
		goCoords:=make([]float64,1)
		return goCoords, fmt.Errorf("Error reading frame")
			}
	goCoords:=make([]float64,totalcoords)
	for j:=0;j<totalcoords;j++{
		goCoords[j]=float64(Ccoords[j])
		}
	return goCoords, nil		
	}

//Similar to XtcGetFrame but takes the buffer given to C as a parameter
//This allows to use the same buffer everytime, simplifying the memory allocation
func xtcGetFrameEfficient(fp *C.XDRFILE, Ccoords []C.float, natoms int)([]float64,error){
	totalcoords:=natoms*3
	cnatoms:=C.int(natoms)
	worked:=C.get_coords(fp,&Ccoords[0],cnatoms)
	if worked==11{
		goCoords:=make([]float64,1)
		return goCoords, fmt.Errorf("No more frames")
		}
	if worked!=0{
		goCoords:=make([]float64,1)
		return goCoords, fmt.Errorf("Error reading frame")
			}
	goCoords:=make([]float64,totalcoords)
	for j:=0;j<totalcoords;j++{
		goCoords[j]=float64(Ccoords[j])
		}
	return goCoords, nil		
	}
	

func XtcGetFrameDrop(fp *C.XDRFILE,natoms int)(error){
	totalcoords:=natoms*3
	cnatoms:=C.int(natoms)
	Ccoords:= make([]C.float,totalcoords)
	worked:=C.get_coords(fp,&Ccoords[0],cnatoms)
	if worked==11{
		return fmt.Errorf("No more frames")
		}
	if worked!=0{
		return fmt.Errorf("Error reading frame")
			}
	return nil		
	}



//Similar to XtcGetFrame but takes the buffer given to C as a parameter
//This allows to use the same buffer everytime, simplifying the memory allocation
//In addition, this doesnt allocate float64 slice, but simply drops the read frame
func xtcGetFrameEfficientDrop(fp *C.XDRFILE, Ccoords []C.float, natoms int)(error){
	cnatoms:=C.int(natoms)
	worked:=C.get_coords(fp,&Ccoords[0],cnatoms)
	if worked==11{
		return fmt.Errorf("No more frames")
		}
	if worked!=0{
		return fmt.Errorf("Error reading frame")
			}
	return nil		
	}
	

func XtcClose(fp *C.XDRFILE){
	C.xtc_close(fp)
	}
