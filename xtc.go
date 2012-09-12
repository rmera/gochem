// +build  gromacs 

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

/*ReadXtcFrames opens the Gromacs trajectory xtc file with name filename
and reads the coordinates for frames starting from ini to end (or the 
last frame in the trajectory) and skipping skip frame between each 
read. The frames are returned as a slice of matrix.Densematrix.
 It returns also the number of frames read, and
 error/nil in failure/success. Note that if there are less frames than
 end, the function wont return error, just the read frames and
 the number of them.*/
func ReadXtcFrames(ini, end, skip int, filename string)([]*matrix.DenseMatrix,int, error) {
	Coords:=make([]*matrix.DenseMatrix,0,1) // I might attempt to give the capacity later
	natoms,err:=XtcCountAtoms(filename)
	if err!=nil{
		return nil, 0, err
		}
	fp,_:=XtcOpen(filename) //We already tested that the name is ok, no need to catch this error
	defer XtcClose(fp)
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
	


//XtcCountAtoms takes the name of a Gromacs xtc trajectory file and returns the 
//number of atoms per frame in the trajectory.
func XtcCountAtoms(name string)(int, error){
	Cfilename:=C.CString(name)
	Cnatoms:=C.read_natoms(Cfilename)
	natoms:=int(Cnatoms)
	return natoms, nil
	}
	
/*XtcOpen Opens the Gromacs xtc trajectory file with the name name and
 * returns a C pointer to it, which can passed to other functions of
 * the package to  */
func  XtcOpen(name string)(*C.XDRFILE, error){
	Cfilename:=C.CString(name)
	fp:=C.openfile(Cfilename)
	if fp==nil{
		return nil, fmt.Errorf("Unable to open xtc file")
		}
	return fp, nil
	}

/*XtcGetFrame takes a C pointer to an open Gromacs xtc file and the 
 * number of atoms per frame in the file. It reads the coordinates
 * for the next frame of the file and returns  them as a slice of
 * float64 */
func XtcGetFrame(fp *C.XDRFILE,natoms int)(*matrix.DenseMatrix,error){
	totalcoords:=natoms*3
	cnatoms:=C.int(natoms)
	Ccoords:= make([]C.float,totalcoords)
	worked:=C.get_coords(fp,&Ccoords[0],cnatoms)
	if worked==11{
		return nil, fmt.Errorf("No more frames") //This is not really an error and should be catched in the calling function
		}
	if worked!=0{
		return nil, fmt.Errorf("Error reading frame")
			}
	goCoords:=make([]float64,totalcoords)
	for j:=0;j<totalcoords;j++{
		goCoords[j]=10*(float64(Ccoords[j]))  //nm to Angstroms
		}
	return matrix.MakeDenseMatrix(goCoords,natoms,3), nil		
	}

/*XtcGetFrameEfficient takes a C pointer to an open Gromacs xtc file, 
 * a slice of C float with enough size to contain all the coordinates
 * to be read, and number of atoms per frame in the file. 
 * It reads the coordinates for the next frame of the file and returns 
 * them as a slice of float64. The fact that it takes the intermediate
 * buffer as an argument means that one can save many memory allocations.
 * It returns nill on success, a "No more frames" error if no more 
 * frames, and other error in other case.*/
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
		goCoords[j]=10*float64(Ccoords[j])  //nm to angstroms
		}
	return goCoords, nil		
	}
	
/*XtcGetFrameDrop takes a C pointer to an open Gromacs xtc file and the 
 * number of atoms per frame in the file. It reads the coordinates
 * for the next frame of the file and discards them, returning only 
 * error/nil in case of failure/success. The special case of no more
 * frames to read causes it to return a "No more frames" error.*/
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


/*XtcGetFrameEfficientDrop takes a C pointer to an open Gromacs xtc file, 
 * a slice of C float with enough size to contain all the coordinates
 * to be read, and number of atoms per frame in the file. 
 * It reads the coordinates for the next frame of the file and discart 
 * them, returning error/nil in case of failure/success. 
 * The fact that it takes the intermediate
 * buffer as an argument means that one can save many memory allocations. */
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
	
/*XtcClose takes a pointer to an open Gromacs xtc trajectory file
 * and closes the file pointed by the pointer.*/
func XtcClose(fp *C.XDRFILE){
	C.xtc_close(fp)
	}
