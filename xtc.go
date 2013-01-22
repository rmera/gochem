// +build xtc

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
import "runtime"

//Container for an GROMACS XTC binary trajectory file.
type XtcObj struct {
	readable bool
	natoms   int
	filename string
	fp       *C.XDRFILE //pointer to the XDRFILE
	goCoords []float64
	cCoords  []C.float
}


func MakeXtc(filename string) (*XtcObj,error) {
	traj:=new(XtcObj)
	if err := traj.initRead(filename);err!=nil{
		return nil,err
		}
	return traj, nil
	
	}


//Returns true if the object is ready to be read from
//false otherwise. IT doesnt guarantee that there is something
//to read.
func (X *XtcObj) Readable() bool {
	if X.readable {
		return true
	}
	return false
}

//InitRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid
func (X *XtcObj) initRead(name string) error {
	Cfilename := C.CString(name)
	Cnatoms := C.read_natoms(Cfilename)
	X.natoms = int(Cnatoms)
	totalcoords := X.natoms * 3
	X.fp = C.openfile(Cfilename)
	if X.fp == nil {
		return fmt.Errorf("Unable to open xtc file")
	}
	//The idea is to reserve less memory, using the same buffer many times.
	X.goCoords = make([]float64, totalcoords, totalcoords)
	X.cCoords = make([]C.float, totalcoords, totalcoords)
	//This should close the file.
	runtime.SetFinalizer(X, func(X *XtcObj) {
		C.xtc_close(X.fp)
	})
	X.readable=true
	return nil
}

//Next Reads the next frame in a XtcObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (X *XtcObj) Next(keep bool) (*matrix.DenseMatrix, error) {
	if !X.Readable() {
		return nil, fmt.Errorf("Traj object uninitialized to read")
	}
	totalcoords := 3 * X.natoms
	cnatoms := C.int(X.natoms)
	worked := C.get_coords(X.fp, &X.cCoords[0], cnatoms)
	if worked == 11 {
		X.readable=false
		return nil, fmt.Errorf("No more frames") //This is not really an error and should be catched in the calling function
	}
	if worked != 0 {
		X.readable=false
		return nil, fmt.Errorf("Error reading frame")
	}
	if keep == true { //collect the trame
		//make sure the buffer is there.
		if X.goCoords == nil {
			X.goCoords = make([]float64, totalcoords, totalcoords)
		}
		for j := 0; j < totalcoords; j++ {
			X.goCoords[j] = 10 * (float64(X.cCoords[j])) //nm to Angstroms
		}
		return matrix.MakeDenseMatrix(X.goCoords, X.natoms, 3), nil
	}
	return nil, nil //Just drop the frame
}

/*NextConc takes a slice of bools and reads as many frames as elements the list has
form the trajectory. The frames are discarted if the corresponding elemetn of the slice
* is false. The function returns a slice of channels through each of each of which 
* a *matrix.DenseMatrix will be transmited*/
func (X *XtcObj) NextConc(frames []bool) ([]chan *matrix.DenseMatrix, error) {
	//For this we dont need this memory
	if X.goCoords != nil {
		X.goCoords = nil
	}
	if X.natoms == 0 {
		return nil, fmt.Errorf("Traj object uninitialized to read")
	}
	framechans := make([]chan *matrix.DenseMatrix, len(frames)) //the slice of chans that will be returned
	totalcoords := X.natoms * 3
	cnatoms := C.int(X.natoms)
	used := false
	for key, val := range frames {
		cCoords := X.cCoords
		/*There seems to be an issue with Go here, although I might of course be wrong
		if I try to used the buffer in X.goCoords, even if I check that its overwritten with the
		* values from this frame, and that is sent to the channel, somewhow the channel gets the
		* a matrix witht he values previously stored in the buffer. The solution here, just not use
		* the buffer, works but wastes the memory of the buffer*/
		goCoords := make([]float64, totalcoords, totalcoords) //X.goCoords
		if used == true {
			cCoords = make([]C.float, totalcoords, totalcoords)
		}
		worked := C.get_coords(X.fp, &cCoords[0], cnatoms)
		//Error handling
		if worked == 11 {
			if used == false {
				X.readable=false
				return nil, fmt.Errorf("No more frames") //This is not really an error and 
			} else { //should be catched in the calling function
				X.readable=false
				return framechans, fmt.Errorf("No more frames") //same
			}
		}
		if worked != 0 {
			X.readable=false
			return nil, fmt.Errorf("Error reading frame")
		}
		//We have to test for used twice to allow allocating for goCoords
		//When the buffer is not going to be used.
		if val == false {
			framechans[key] = nil //ignored frame
			continue
		}
		//		if used==true{
		//			goCoords=make([]float64,totalcoords,totalcoords) 
		//			}
		used = true
		framechans[key] = make(chan *matrix.DenseMatrix)
		//Now the parallel part
		go func(natoms int, cCoords []C.float, goCoords []float64, pipe chan *matrix.DenseMatrix) {
			for j := 0; j < natoms*3; j++ {
				//				fmt.Println(10*(float64(cCoords[j])))
				goCoords[j] = 10 * (float64(cCoords[j])) //nm to Angstroms
			}
			temp := matrix.MakeDenseMatrix(goCoords, natoms, 3)
			//			fmt.Println("in gorutine!", temp.GetRowVector(2))
			pipe <- temp
		}(X.natoms, cCoords, goCoords, framechans[key])
	}
	return framechans, nil
}


//Natoms returns the number of atoms per frame in the XtcObj.
//XtcObj must be initialized. 0 means an uninitialized object.
func (X *XtcObj) Len() int {
	return X.natoms
}


