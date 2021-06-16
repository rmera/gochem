/*
 * crd.go, part of gochem
 *
 * Copyright 2018 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
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

package amberold

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

//Container for an old-Amber/pDynamo trajectory file.
type CrdObj struct {
	natoms    int
	readLast  bool //Have we read the last frame?
	readable  bool //Is it ready to be read?
	filename  string
	nnew      bool     //Still no frame read from it?
	fixed     int32    //Fixed atoms (not supported)
	ioread    *os.File //The crd file
	crd       *bufio.Reader
	remaining []float64
	box       bool
}

//New creates a new Old Amber trajectory object from a file.
func New(filename string, ats int, box bool) (*CrdObj, error) {
	var err error
	traj := new(CrdObj)
	traj.ioread, err = os.Open(filename)
	if err != nil {
		return nil, Error{UnableToOpen, filename, []string{"New"}, true}

	}
	traj.filename = filename
	traj.crd = bufio.NewReader(traj.ioread)
	_, err = traj.crd.ReadString('\n') //The first line is just a comment
	//	println(TEST) ///////////
	//	TEST,err=traj.crd.ReadString('\n') //The first line is just a comment
	//	println(TEST) ////////////
	if err != nil {
		return nil, err //CHANGE
	}
	traj.natoms = ats
	traj.remaining = make([]float64, 0, 9)
	if box {
		traj.box = true
	}
	traj.readable = true
	return traj, nil
}

//Readable returns true if the object is ready to be read from
//false otherwise. It doesnt guarantee that there is something
//to read.
func (C *CrdObj) Readable() bool {
	return C.readable
}

//Next Reads the next frame in a DcDObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (C *CrdObj) Next(keep *v3.Matrix) error {
	const ncoords = 3
	if !C.readable {
		return Error{TrajUnIni, C.filename, []string{"Next"}, true}
	}
	//What do we do with remaining?
	cont := -1
	line := 0
	//The following allows discarding a frame while still keeping track of it.
	//Everything is the same if you read or discard, except the function that
	//would set the values to the matrix simply does nothing in the discard case.
	var setter func(file, col int, val float64)
	if keep != nil {
		setter = func(file, col int, val float64) {
			//			println(keep.NVecs()) ///////////////////////////////////////////
			keep.Set(file, col, val)

		}
	} else {
		setter = func(file, col int, val float64) {}
	}
	//Here we read the coords that were left remaining in the last read.
	for _, coord := range C.remaining {
		if cont < ncoords-1 {
			cont++
			setter(line, cont, coord)
			continue
		}
		if cont >= ncoords-1 && line < C.natoms-1 {
			line++
			cont = 0
			setter(line, cont, coord)
			continue
		}
	}
	cont = -1
	C.remaining = C.remaining[0:0] //This might not work as expected. I need it to set C.remaining to zero length.
	for line < C.natoms-1 || cont < 2 {
		i, err := C.crd.ReadString('\n')
		//here we assume the error is an EOF. I need to change this to actually check.
		if err != nil {
			C.readable = false
			//			println(err.Error()) //////
			return newlastFrameError(C.filename, "Next")
		}
		l := strings.Fields(i)
		for _, j := range l {
			coord, err := strconv.ParseFloat(j, 64)
			if err != nil {
				return Error{fmt.Sprint("Unable to read coordinates from Amber trajectory", err.Error()), C.filename, []string{"strconv.ParseFloat", "Next"}, true}
			}
			if cont < ncoords-1 {
				cont++
				setter(line, cont, coord)
				continue
			}
			if cont >= ncoords-1 && line < C.natoms-1 {
				line++
				cont = 0
				setter(line, cont, coord)
				continue
			}
			if cont >= cont-1 && line >= C.natoms-1 {

				//				println("ql")////
				C.remaining = append(C.remaining, coord)
			}
		}
		//	println("wei") ////////////
	}
	if C.box {
		C.nextBox()
	}
	return nil

}

func (C *CrdObj) nextVelBox() error {
	var keep *v3.Matrix //nil
	const ncoords = 3
	if !C.readable {
		return Error{TrajUnIni, C.filename, []string{"Next"}, true}
	}
	var natoms int = C.natoms + 1
	//What do we do with remaining?
	cont := 0
	line := 0
	//The following allows discarding a frame while still keeping track of it.
	//Everything is the same if you read or discard, except the function that
	//would set the values to the matrix simply does nothing in the discard case.
	var setter func(file, col int, val float64)
	if keep != nil {
		setter = func(file, col int, val float64) { keep.Set(file, col, val) }
	} else {
		setter = func(file, col int, val float64) {}
	}
	//Here we read the coords that were left remaining in the last read.
	for _, coord := range C.remaining {
		if cont < ncoords {
			setter(line, cont, coord)
			coord++
			continue
		}
		if cont >= ncoords && line < natoms-1 {
			line++
			cont = 0
			setter(line, cont, coord)
			cont++
			continue
		}
	}
	C.remaining = C.remaining[0:0] //This might not work as expected. I need it to set C.remaining to zero length.
	//	println("remaining:", len(C.remaining))
	for line < natoms-1 {
		i, err := C.crd.ReadString('\n')
		//here we assume the error is an EOF. I need to change this to actually check.
		if err != nil {
			C.readable = false
			//			println(err.Error()) //////
			return newlastFrameError(C.filename, "Next")
		}
		l := strings.Fields(i)
		//	fmt.Println(l) ////////////////////////////////////////////////////////////
		//	println("no") ///////////////////////
		const ncoords = 3
		for _, j := range l {
			coord, err := strconv.ParseFloat(j, 64)
			if err != nil {
				return Error{fmt.Sprint("Unable to read coordinates from Amber trajectory", err.Error()), C.filename, []string{"strconv.ParseFloat", "Next"}, true}
			}
			if cont < ncoords {
				//	println("no")////
				setter(line, cont, coord)
				cont++
				continue
			}
			if cont >= ncoords && line < natoms-1 {
				//		println("wei")////
				line++
				cont = 0
				setter(line, cont, coord)
				cont++
				continue
			}
			if cont >= cont && line >= natoms-1 {

				//				println("ql")////
				C.remaining = append(C.remaining, coord)
			}
		}
		//	println("wei") ////////////
	}
	return nil

}

func (C *CrdObj) nextBox() error {
	var keep *v3.Matrix //nil
	const ncoords = 3
	if !C.readable {
		return Error{TrajUnIni, C.filename, []string{"Next"}, true}
	}
	var natoms int = 1
	//What do we do with remaining?
	cont := 0
	line := 0
	//The following allows discarding a frame while still keeping track of it.
	//Everything is the same if you read or discard, except the function that
	//would set the values to the matrix simply does nothing in the discard case.
	var setter func(file, col int, val float64)
	if keep != nil {
		setter = func(file, col int, val float64) { keep.Set(file, col, val) }
	} else {
		setter = func(file, col int, val float64) {}
	}
	//Here we read the coords that were left remaining in the last read.
	for _, coord := range C.remaining {
		if cont < ncoords {
			setter(line, cont, coord)
			coord++
			continue
		}
		if cont >= ncoords && line < natoms-1 {
			line++
			cont = 0
			setter(line, cont, coord)
			cont++
			continue
		}
	}
	C.remaining = C.remaining[0:0] //This might not work as expected. I need it to set C.remaining to zero length.
	//	println("remaining:", len(C.remaining))
	for line < natoms-1 {
		i, err := C.crd.ReadString('\n')
		//here we assume the error is an EOF. I need to change this to actually check.
		if err != nil {
			C.readable = false
			//			println(err.Error()) //////
			return newlastFrameError(C.filename, "Next")
		}
		l := strings.Fields(i)
		//	fmt.Println(l) ////////////////////////////////////////////////////////////
		//	println("no") ///////////////////////
		const ncoords = 3
		for _, j := range l {
			coord, err := strconv.ParseFloat(j, 64)
			if err != nil {
				return Error{fmt.Sprint("Unable to read coordinates from Amber trajectory", err.Error()), C.filename, []string{"strconv.ParseFloat", "Next"}, true}
			}
			if cont < ncoords {
				//	println("no")////
				setter(line, cont, coord)
				cont++
				continue
			}
			if cont >= ncoords && line < natoms-1 {
				//		println("wei")////
				line++
				cont = 0
				setter(line, cont, coord)
				cont++
				continue
			}
			if cont >= cont && line >= natoms-1 {

				//				println("ql")////
				C.remaining = append(C.remaining, coord)
			}
		}
		//	println("wei") ////////////
	}
	return nil

}

//Natoms returns the number of atoms per frame in the XtcObj.
//XtcObj must be initialized. 0 means an uninitialized object.
func (C *CrdObj) Len() int {
	return int(C.natoms)
}

//Errors

//errDecorate is a helper function that asserts that the error is
//implements chem.Error and decorates the error with the caller's name before returning it.
//if used with a non-chem.Error error, it will cause a panic.
func errDecorate(err error, caller string) error {
	err2 := err.(chem.Error) //I know that is the type returned byt initRead
	err2.Decorate(caller)
	return err2
}

//Error is the general structure for Crd trajectory errors. It fullfills  chem.Error and chem.TrajError
type Error struct {
	message  string
	filename string //the input file that has problems, or empty string if none.
	deco     []string
	critical bool
}

func (err Error) Error() string {
	return fmt.Sprintf("Old Amber trajectory file %s error: %s", err.filename, err.message)
}

func (E Error) Decorate(deco string) []string {
	//Even thought this method does not use a pointer as a receiver, and tries to alter the received,
	//it should work, since E.deco is a slice, and hence a pointer itself.

	if deco != "" {
		E.deco = append(E.deco, deco)
	}
	return E.deco
}

func (err Error) FileName() string { return err.filename }

func (err Error) Format() string { return "Old Amber" }

func (err Error) Critical() bool { return err.critical }

const (
	TrajUnIni           = "Traj object uninitialized to read"
	ReadError           = "Error reading frame"
	UnableToOpen        = "Unable to open file"
	SecurityCheckFailed = "FailedSecurityCheck"
	WrongFormat         = "Wrong format in the trajectory file or frame"
	NotEnoughSpace      = "Not enough space in passed blocks"
	EOF                 = "EOF"
)

//lastFrameError implements chem.LastFrameError
type lastFrameError struct {
	deco     []string
	fileName string
}

//lastFrameError does nothing
func (E lastFrameError) NormalLastFrameTermination() {}

func (E lastFrameError) FileName() string { return E.fileName }

func (E lastFrameError) Error() string { return "EOF" }

func (E lastFrameError) Critical() bool { return false }

func (E lastFrameError) Format() string { return "Old Amber" }

func (E lastFrameError) Decorate(deco string) []string {
	//Even thought this method does not use a pointer as a receiver, and tries to alter the received,
	//it should work, since E.deco is a slice, and hence a pointer itself.
	if deco != "" {
		E.deco = append(E.deco, deco)
	}
	return E.deco
}

func newlastFrameError(filename string, caller string) *lastFrameError {
	e := new(lastFrameError)
	e.fileName = filename
	e.deco = []string{caller}
	return e
}
