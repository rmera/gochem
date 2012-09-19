// +build  gromacs 

/*
 * xtc.go, part of gochem
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
//import "unsafe"
import "runtime"


/*The plan is equate PDBs XTCs and in the future DCDs. One needs to separate the molecule methods between actual molecule methods, that requires atoms and coordinates, []atom methods, and  DenseMatrix
 * methods. Then one must implements objects for Xtc trajs that are more or less equivalent to molecules and set up an interface so many analyses can be carried out exactly the same from
 * multiPDB or XTC or (eventually) DCD*/

//Traj is an interface for any trajectory object, including a Molecule Object
type Traj interface{
	//Opens the file and prepares for reading
	InitRead(string) error
	//reads the next frame and returns it as DenseMatrix if keep==true, or discards it if false
	Next(keep bool) *matrix.DenseMatrix
	//Read frames from Traj from ini to end skipping skip frames between read. Returns a slice with coords of each frame
	//the number of frames read and error or nil.
	ManyFrames(ini, end, skip int)([]*matrix.DenseMatrix,int, error)
	}

//Container for an GROMACS XTC binary trajectory file.
type XtcObj struct{
	natoms int
	filename string
	fp *C.XDRFILE //pointer to the XDRFILE
	goCoords []float64
	cCoords []C.float
	}

//InitRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid
func (X *XtcObj) InitRead(name string) error{
	Cfilename:=C.CString(name)
	Cnatoms:=C.read_natoms(Cfilename)
	X.natoms=int(Cnatoms)
	totalcoords:=X.natoms*3
	X.fp=C.openfile(Cfilename)
	if X.fp==nil{
		return fmt.Errorf("Unable to open xtc file")
		}
	//The idea is to reserve less memory, using the same buffer many times.
	X.goCoords=make([]float64,totalcoords,totalcoords)
	X.cCoords= make([]C.float,totalcoords,totalcoords)
	//This should close the file.
	runtime.SetFinalizer(X,func (X *XtcObj){
		C.xtc_close(X.fp)
		})
	return nil
	}

//Next Reads the next frame in a XtcObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (X *XtcObj)Next(keep bool)(*matrix.DenseMatrix, error){
	if X.natoms==0{
		return nil, fmt.Errorf("Traj object uninitialized!")
		}
	totalcoords:=3*X.natoms
	cnatoms:=C.int(X.natoms)
	worked:=C.get_coords(X.fp,&X.cCoords[0],cnatoms)
	if worked==11{
		return nil, fmt.Errorf("No more frames") //This is not really an error and should be catched in the calling function
		}
	if worked!=0{
		return nil, fmt.Errorf("Error reading frame")
			}
	if keep==true{ //collect the trame
		for j:=0;j<totalcoords;j++{
			X.goCoords[j]=10*(float64(X.cCoords[j]))  //nm to Angstroms
			}
		return matrix.MakeDenseMatrix(X.goCoords,X.natoms,3), nil		
		}
	return nil, nil //Just drop the frame
	}

/*NextConc takes a slice of bools and reads as many frames as elements the list has
form the trajectory. The frames are discarted if the corresponding elemetn of the slice
* is false. The function returns a slice of channels through each of each of which 
* a *matrix.DenseMatrix will be transmited*/
func (X *XtcObj)NextConc(frames []bool)([]chan *matrix.DenseMatrix, error){
	//this function is rather ugly since its tuned to performance.
	framechans:=make([]chan *matrix.DenseMatrix,0,len(frames))  //the slice of chans that will be returned
	used:=false//whether we have actually read a frame (and not dropped them)
	totalcoords:=X.natoms*3
	cnatoms:=C.int(X.natoms)
	for _,val:=range(frames){  //We use the buffers in the traj object when possible.
		cCoords:=X.cCoords   
		goCoords:=X.goCoords
		if used==true{  //If we have the previous buffers in use we need to allocate new memory.
			cCoords=make([]C.float,totalcoords,totalcoords) //the memory for goCoords is allocated later if the frame is not dropped, see MARK.
			}
		worked:=C.get_coords(X.fp,&X.cCoords[0],cnatoms)
		if val==false{
			framechans=append(framechans,nil)  //ignored frame
			continue
			}
		framechans=append(framechans,make(chan *matrix.DenseMatrix))
		if used==true{
			goCoords=make([]float64,totalcoords,totalcoords)   //MARK
			}
		used=true //Marks that we have already read at least one frame.
		//Error handling
		if worked==11{
			return nil, fmt.Errorf("No more frames") //This is not really an error and should be catched in the calling function
			}
		if worked!=0{
			return nil, fmt.Errorf("Error reading frame")
				}
		//Now the parallel part
		go func(natoms int,cCoords []C.float, goCoords []float64, pipe chan *matrix.DenseMatrix){
			for j:=0;j<natoms*3;j++{
				goCoords[j]=10*(float64(cCoords[j]))  //nm to Angstroms
				}
			pipe<-matrix.MakeDenseMatrix(goCoords,natoms,3)				
			}(X.natoms,cCoords,goCoords,framechans[len(framechans)-1])
		}
	return framechans,nil
	}


//Read frames from Traj from ini to end skipping skip frames between read. Returns a slice with coords of each frame
//the number of frames read and error or nil.
func (X *XtcObj)ManyFrames(ini, end, skip int)([]*matrix.DenseMatrix,int, error) {
	Coords:=make([]*matrix.DenseMatrix,0,1) // I might attempt to give the capacity later
	i:=0
	for ;;i++{
		if i>end {
			break
			}
		if i<ini || (i-(1+ini))%skip!=0{
			_,err:=X.Next(false) //Drop the frame
			if err!=nil{
				if err.Error()=="No more frames"{
					break //No more frames is not really an error
					}else{
					return Coords, i, err	
					}
				}
			}else{
			coords,err:=X.Next(true)
			if err!=nil{
				if err.Error()=="No more frames"{
					break //No more frames is not really an error
					}else{
					return Coords, i, err
					}
				}
			Coords=append(Coords,coords)
			}
		}
	return Coords, i, nil
	}








	

