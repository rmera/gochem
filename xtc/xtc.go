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

package xtc

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
import "runtime"
import "github.com/rmera/gochem"

//Container for an GROMACS XTC binary trajectory file.
type XTCObj struct {
	readable   bool
	natoms     int
	filename   string
	fp         *C.XDRFILE //pointer to the XDRFILE
	goCoords   []float64
	concBuffer [][]C.float
	cCoords    []C.float
	buffSize   int
}

func New(filename string) (*XTCObj, error) {
	traj := new(XTCObj)
	if err := traj.initRead(filename); err != nil {
		return nil, err
	}
	return traj, nil

}

//Returns true if the object is ready to be read from
//false otherwise. IT doesnt guarantee that there is something
//to read.
func (X *XTCObj) Readable() bool {
	if X.readable {
		return true
	}
	return false
}

//InitRead initializes a XTCObj for reading.
//It requires only the filename, which must be valid
func (X *XTCObj) initRead(name string) error {
	Cfilename := C.CString(name)
	Cnatoms := C.read_natoms(Cfilename)
	X.natoms = int(Cnatoms)
	totalcoords := X.natoms * 3
	X.fp = C.openfile(Cfilename)
	if X.fp == nil {
		return fmt.Errorf("Unable to open xtc file")
	}
	//The idea is to reserve less memory, using the same buffer many times.
	//X.goCoords = make([]float64, totalcoords, totalcoords)
	X.cCoords = make([]C.float, totalcoords, totalcoords)
	X.concBuffer = append(X.concBuffer, X.cCoords)
	X.buffSize = 1
	//This should close the file.
	runtime.SetFinalizer(X, func(X *XTCObj) {
		C.xtc_close(X.fp)
	})
	X.readable = true
	return nil
}

//Next Reads the next frame in a XTCObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (X *XTCObj) Next(output *chem.VecMatrix) error {
	if !X.Readable() {
		return fmt.Errorf("Traj object uninitialized to read")
	}
	cnatoms := C.int(X.natoms)
	worked := C.get_coords(X.fp, &X.cCoords[0], cnatoms)
	if worked == 11 {
		X.readable = false
		return newlastFrameError(X.filename) //This is not really an error and should be catched in the calling function
	}
	if worked != 0 {
		X.readable = false
		return fmt.Errorf("Error reading frame")
	}
	if output != nil { //col the frame
		r, c := output.Dims()
		if r < (X.natoms) {
			panic("VecMatrix too small!")
		}
		for j := 0; j < r; j++ {
			for k := 0; k < c; k++ {
				l := k + (3 * j)
				output.Set(j, k, (10 * float64(X.cCoords[l]))) //nm to Angstroms
			}
		}
		return nil
	}
	return nil //Just drop the frame
}

//SetConcBuffer
func (X *XTCObj) setConcBuffer(batchsize int) error {
	l := X.buffSize
	if l == batchsize {
		return nil
	} else if l > batchsize {
		for i := batchsize; i < l; i++ {
			X.concBuffer[i] = nil //no idea if this actually works
		}
		X.concBuffer = X.concBuffer[:batchsize-1] //not sure if this frees the remaining []float32 slices
		X.buffSize = batchsize
		return nil
	}
	for i := 0; i < batchsize-l; i++ {
		tmp := make([]C.float, X.Len()*3)
		X.concBuffer = append(X.concBuffer, tmp)
	}
	X.buffSize = batchsize
	return nil
}

/*NextConc takes a slice of bools and reads as many frames as elements the list has
form the trajectory. The frames are discarted if the corresponding elemetn of the slice
* is false. The function returns a slice of channels through each of each of which
* a *matrix.DenseMatrix will be transmited*/
func (X *XTCObj) NextConc(frames []*chem.VecMatrix) ([]chan *chem.VecMatrix, error) {
	if X.buffSize < len(frames) {
		X.setConcBuffer(len(frames))
	}
	if X.natoms == 0 {
		return nil, fmt.Errorf("Traj object uninitialized to read")
	}
	framechans := make([]chan *chem.VecMatrix, len(frames)) //the slice of chans that will be returned
	cnatoms := C.int(X.natoms)
	used := false
	for key, val := range frames {
		//cCoords:=X.concBuffer[key]
		worked := C.get_coords(X.fp, &X.concBuffer[key][0], cnatoms)
		//Error handling
		if worked == 11 {
			if used == false {
				X.readable = false
				return nil, newlastFrameError(X.filename) //This is not really an error and
			} else { //should be catched in the calling function
				X.readable = false
				return framechans, newlastFrameError(X.filename) //same
			}
		}
		if worked != 0 {
			X.readable = false
			return nil, fmt.Errorf("Error reading frame")
		}
		if val == nil {
			framechans[key] = nil //ignored frame
			continue
		}
		used = true
		framechans[key] = make(chan *chem.VecMatrix)
		//Now the parallel part
		go func(natoms int, cCoords []C.float, goCoords *chem.VecMatrix, pipe chan *chem.VecMatrix) {
			r, c := goCoords.Dims()
			for j := 0; j < r; j++ {
				for k := 0; k < c; k++ {
					l := k + (3 * j)
					goCoords.Set(j, k, (10 * float64(cCoords[l]))) //nm to Angstroms
				}
			}
			pipe <- goCoords
		}(X.natoms, X.concBuffer[key], val, framechans[key])
	}
	return framechans, nil
}

//Natoms returns the number of atoms per frame in the XTCObj.
//XTCObj must be initialized. 0 means an uninitialized object.
func (X *XTCObj) Len() int {
	return X.natoms
}


//Errors

type lastFrameError struct {
	fileName string
}


func (E *lastFrameError) Error() string {
	return fmt.Sprintf("No More Frames: Last frame in xtc trajectory from file %10s reached", E.fileName)
}

func (E *lastFrameError) Format() string {
	return "xtc"
}


func (E *lastFrameError) FileName() string{
	return E.fileName
}

func (E *lastFrameError) NormalLastFrameTermination(){
}



func newlastFrameError(filename string) *lastFrameError{
	e:=new(lastFrameError)
	e.fileName=filename
	return e
}

