/*
 * dcd.go, part of gochem
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

package dcd

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"os"
	"runtime"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

const mAXTITLE int32 = 80
const rSCAL32BITS int32 = 1

//Container for an Charmm/NAMD binary trajectory file.
type DCDObj struct {
	natoms     int32
	buffSize   int
	readLast   bool //Have we read the last frame?
	readable   bool //Is it ready to be read?
	filename   string
	charmm     bool //Charmm traj?
	extrablock bool
	fourdim    bool
	new        bool     //Still no frame read from it?
	fixed      int32    //Fixed atoms (not supported)
	dcd        *os.File //The DCD file
	dcdFields  [][]float32
	concBuffer [][][]float32
	endian     binary.ByteOrder
}

//New builds a new DCDObj object from a DCD trajectory file
func New(filename string) (*DCDObj, error) {
	traj := new(DCDObj)
	if err := traj.initRead(filename); err != nil {
		return nil, errDecorate(err, "New")
	}
	traj.dcdFields = make([][]float32, 3, 3)
	traj.dcdFields[0] = make([]float32, int(traj.natoms), int(traj.natoms))
	traj.dcdFields[1] = make([]float32, int(traj.natoms), int(traj.natoms))
	traj.dcdFields[2] = make([]float32, int(traj.natoms), int(traj.natoms))
	traj.concBuffer = append(traj.concBuffer, traj.dcdFields)
	return traj, nil

}

//Readable returns true if the object is ready to be read from
//false otherwise. It doesnt guarantee that there is something
//to read.
func (D *DCDObj) Readable() bool {
	return D.readable
}

//initRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid.
//It support big and little endianness, charmm or (namd>=2.1) and no
//fixed atoms.
func (D *DCDObj) initRead(name string) error {
	wrapbinerr := func(err error) error {
		return Error{err.Error(), D.filename, []string{"binary.Read", "initRead"}, true}
	}

	rec_scale := rSCAL32BITS //At least for now we will not support anything else.
	D.endian = binary.LittleEndian
	_ = rec_scale
	NB := bytes.NewBuffer //shortness sake
	var err error
	D.dcd, err = os.Open(name)
	if err != nil {
		return Error{err.Error(), D.filename, []string{"os.Open", "initRead"}, true}
	}
	var check int32
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return wrapbinerr(err)
	}
	//For some reason the first thing we should read is an 84.
	//If this fails it means that the file is big endian.
	if check != 84 {
		D.endian = binary.BigEndian //
	}
	//Then the magic number "CORD", also for some unknown reason.
	magic := make([]byte, 4, 4)
	if err := binary.Read(D.dcd, D.endian, magic); err != nil {
		return wrapbinerr(err)
	}
	if string(magic) != "CORD" {
		return Error{WrongFormat + ": Wrong magic number", D.filename, []string{"initRead"}, true}
	}

	//We first read a big chuck for random access.
	buf := make([]byte, 80, 80)
	if err := binary.Read(D.dcd, D.endian, buf); err != nil {
		return wrapbinerr(err)
	}
	//X-plor sets this last int to zero, charmm sets it to its version number.
	//if we have a charmm file we get some additional flags.
	if err := binary.Read(NB(buf[76:]), D.endian, &check); err != nil {
		return wrapbinerr(err)
	}
	if check != 0 {
		//		fmt.Println("CHARMM!!!") //////77
		D.charmm = true
		if err := binary.Read(NB(buf[40:]), D.endian, &check); err != nil {
			return wrapbinerr(err)
		}
		if check != 0 {
			//			fmt.Println("block", check) ///////////
			D.extrablock = true
		}
		if err := binary.Read(NB(buf[40:]), D.endian, &check); err != nil {
			return wrapbinerr(err)
		}
		if check == 1 {
			//			fmt.Println("4-dim", check) ///////////
			D.fourdim = true
		}

	} else {
		return Error{WrongFormat + ": X-plor DCD not supported", D.filename, []string{"initRead"}, true}
	}
	if err := binary.Read(NB(buf[32:]), D.endian, &D.fixed); err != nil {
		return Error{err.Error(), D.filename, []string{"initRead"}, true}

	}
	//	fmt.Println("fixed", D.fixed)
	var delta float32 //This should work only on Charmm and namd >=2.1
	if err := binary.Read(NB(buf[36:]), D.endian, &delta); err != nil {
		return wrapbinerr(err)
	}
	//	fmt.Println("delta:", delta)///////////////////////////////////////

	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return Error{err.Error(), D.filename, []string{"initRead"}, true}
	}
	if check != 84 {
		return Error{WrongFormat, D.filename, []string{"initRead"}, true}
	}
	var input_int int32
	if err := binary.Read(D.dcd, D.endian, &input_int); err != nil {
		return wrapbinerr(err)
	}
	//how many units of mAXTITLE does the title have?
	var ntitle int32
	if err := binary.Read(D.dcd, D.endian, &ntitle); err != nil {
		return wrapbinerr(err)
	}
	title := make([]byte, mAXTITLE*ntitle, mAXTITLE*ntitle)
	if err := binary.Read(D.dcd, D.endian, title); err != nil {
		return wrapbinerr(err)
	}
	//	fmt.Println("Title:", string(title))///////////////////////////////////////
	if err := binary.Read(D.dcd, D.endian, &input_int); err != nil {
		return wrapbinerr(err)

	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return wrapbinerr(err)

	}
	if check != 4 { //one must read a 4 before the natoms
		return Error{WrongFormat, D.filename, []string{"initRead"}, true}
	}
	if err := binary.Read(D.dcd, D.endian, &D.natoms); err != nil {
		return wrapbinerr(err)
	}
	//	fmt.Println("natoms", D.natoms)
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return wrapbinerr(err)
	}
	if check != 4 { //and one more 4
		return Error{WrongFormat, D.filename, []string{"initRead"}, true}
	}
	if D.fixed == 0 {
		runtime.SetFinalizer(D, func(D *DCDObj) {
			D.dcd.Close()
		})
		D.readable = true
		return nil //nothing else to do
	}
	D.new = true //nothing read yet
	return Error{"Fixed atoms not supported", D.filename, []string{"initRead"}, true}

}

//Next Reads the next frame in a DcDObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DCDObj) Next(keep *v3.Matrix) error {
	if !D.readable {
		return Error{TrajUnIni, D.filename, []string{"Next"}, true}
	}
	if D.dcdFields == nil {
		D.dcdFields = make([][]float32, 3, 3)
		D.dcdFields[0] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[1] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[2] = make([]float32, int(D.natoms), int(D.natoms))
	}
	if err := D.nextRaw(D.dcdFields); err != nil {
		return errDecorate(err, "Next")
	}
	if keep == nil {
		return nil
	}
	if r, _ := keep.Dims(); int32(r) < D.natoms {
		panic("Not enough space in matrix")
	}
	//outBlock := make([]float64, int(D.natoms)*3, int(D.natoms)*3)
	for i := 0; i < int(D.natoms); i++ {
		k := i - i*(i/int(D.natoms))
		keep.Set(k, 0, float64(D.dcdFields[0][k]))
		keep.Set(k, 1, float64(D.dcdFields[1][k]))
		keep.Set(k, 2, float64(D.dcdFields[2][k]))
	}
	//	final := NewVecs(outBlock, int(D.natoms), 3)
	//	fmt.Print(final)/////////7
	return nil
}

//Next Reads the next frame in a XtcObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DCDObj) nextRaw(blocks [][]float32) error {
	if len(blocks[0]) != int(D.natoms) || len(blocks[1]) != int(D.natoms) || len(blocks[2]) != int(D.natoms) {
		return Error{NotEnoughSpace, D.filename, []string{"nextRaw"}, true}
	}
	D.new = false
	if D.readLast {
		D.readable = false
		return newlastFrameError(D.filename, "nextRaw")
	}
	//if there is an extra block we just skip it.
	//Sadly, even when there is an extra block, it is not present in all
	//snapshots for some trajectories, so we must use the block size to see if
	//there is an extra block or if the X block starts inmediately
	var blocksize int32
	if D.extrablock {
		if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
			return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}
		}
		//If the blocksize is 4*natoms it means that the block is not an
		//extra block, but the X coordinates, and thus we must skip the following
		if blocksize != D.natoms*4 {
			if _, err := D.readByteBlock(blocksize); err != nil {
				return err
			}
			blocksize = 0
		}
	}
	//now get the coords, each as a slice of float32
	//X
	//we collect the X block size again only if it has not been collected before
	if blocksize == 0 {
		if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
			return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}

		}
	}
	err := D.readFloat32Block(blocksize, blocks[0])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//	fmt.Println("X", len(xblock)) //, xblock)
	//	fmt.Println("X", blocks[0])  ///////////////////////////////
	//Y
	//Collect the size first, then the rest
	if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
		return err
	}
	err = D.readFloat32Block(blocksize, blocks[1])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//	fmt.Println("Y", blocks[1])
	//Z
	if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}
	}
	err = D.readFloat32Block(blocksize, blocks[2])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//	fmt.Println("Z", blocks[2])
	//we skip the 4-D values if they exist. Apparently this is not present in the
	//last snapshot, so we use an EOF here to signal that we have read the last snapshot.
	if D.charmm && D.fourdim {
		if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
			if err.Error() == "EOF" {
				D.readLast = true
				//	fmt.Println("LAST!")
			} else {
				return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}
			}
		}
		if !D.readLast {
			if _, err := D.readByteBlock(blocksize); err != nil {
				return errDecorate(err, "nextRaw")
			}
		}
	}
	return nil

}

//Queries the size of a block, and reads its contents into block, which must have the
//appropiate size.
func (D *DCDObj) readFloat32Block(blocksize int32, block []float32) error {
	var check int32
	//	fmt.Println("blockf",blocksize)
	if err := binary.Read(D.dcd, D.endian, block); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Read", "readFloat32Block"}, true}

	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Read", "readFloat32Block"}, true}
	}
	if check != blocksize {
		return Error{WrongFormat, D.filename, []string{"readFloat32Block"}, true}
	}
	return nil
}

//Queries the size of a block, make a slice of a quarter of that size
//and reads that ammount of float32. This function is used for the
//
func (D *DCDObj) readByteBlock(blocksize int32) ([]byte, error) {
	var check int32
	//	fmt.Println("blockb",blocksize)
	block := make([]byte, blocksize, blocksize)
	if err := binary.Read(D.dcd, D.endian, block); err != nil {
		return nil, Error{err.Error(), D.filename, []string{"binary.Read", "readByteBlock"}, true}

	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return nil, Error{err.Error(), D.filename, []string{"binary.Read", "readByteBlock"}, true}

	}
	if check != blocksize {
		return nil, Error{SecurityCheckFailed, D.filename, []string{"readByteBlock"}, true}
	}
	return block, nil
}

//Len returns the number of atoms per frame in the XtcObj.
//XtcObj must be initialized. 0 means an uninitialized object.
func (D *DCDObj) Len() int {
	return int(D.natoms)
}

//This function never actually returns error. Still, it is an internal function
//and I may modify it some days so I won't change the signature.
func (D *DCDObj) setConcBuffer(batchsize int) error {
	l := D.buffSize
	if l == batchsize {
		return nil
	} else if l > batchsize {
		for i := batchsize; i < l; i++ {
			for j, _ := range D.concBuffer[i] {
				D.concBuffer[i][j] = nil
			}
			D.concBuffer[i] = nil //no idea if this actually works. Edit: According to tests it does work.
		}
		D.concBuffer = D.concBuffer[:batchsize-1] //not sure if this frees the remaining []float32 slices
		D.buffSize = batchsize
		return nil
	}
	for i := 0; i < batchsize-l; i++ {
		x := make([]float32, D.Len())
		y := make([]float32, D.Len())
		z := make([]float32, D.Len())
		tmp := [][]float32{x, y, z}
		D.concBuffer = append(D.concBuffer, tmp)
	}
	D.buffSize = batchsize
	return nil
}

//NextConc takes a slice of bools and reads as many frames as elements the list has
//form the trajectory. The frames are discarted if the corresponding elemetn of the slice
//is false. The function returns a slice of channels through each of each of which
// a *matrix.DenseMatrix will be transmited
func (D *DCDObj) NextConc(frames []*v3.Matrix) ([]chan *v3.Matrix, error) {
	if !D.Readable() {
		return nil, Error{TrajUnIni, D.filename, []string{"NextConc"}, true}
	}
	framechans := make([]chan *v3.Matrix, len(frames)) //the slice of chans that will be returned
	if D.buffSize < len(frames) {
		D.setConcBuffer(len(frames))
	}
	for key, _ := range frames {
		DFields := D.concBuffer[key]
		if err := D.nextRaw(DFields); err != nil {
			return nil, errDecorate(err, "NextConc")
		}
		//We have to test for used twice to allow allocating for goCoords
		//When the buffer is not going to be used.
		if frames[key] == nil {
			framechans[key] = nil //ignored frame
			continue
		}
		framechans[key] = make(chan *v3.Matrix)
		//Now the parallel part
		go func(natoms int, DFields [][]float32, keep *v3.Matrix, pipe chan *v3.Matrix) {
			for i := 0; i < int(D.natoms); i++ {
				k := i - i*(i/int(D.natoms))
				keep.Set(k, 0, float64(DFields[0][k]))
				keep.Set(k, 1, float64(DFields[1][k]))
				keep.Set(k, 2, float64(DFields[2][k]))
			}
			//			fmt.Println("in gorutine!", temp.GetRowVector(2))
			pipe <- keep
		}(int(D.natoms), DFields, frames[key], framechans[key])
	}
	return framechans, nil
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

//Error is the general structure for DCD trajectory errors. It fullfills  chem.Error and chem.TrajError
type Error struct {
	message  string
	filename string //the input file that has problems, or empty string if none.
	deco     []string
	critical bool
}

func (err Error) Error() string {
	return fmt.Sprintf("dcd file %s error: %s", err.filename, err.message)
}

//Decorate Adds new information to the error
func (E Error) Decorate(deco string) []string {
	//Even thought this method does not use a pointer as a receiver, and tries to alter the received,
	//it should work, since E.deco is a slice, and hence a pointer itself.

	if deco != "" {
		E.deco = append(E.deco, deco)
	}
	return E.deco
}

//Filename returns the file to which the failing trajectory was associated
func (err Error) FileName() string { return err.filename }

//Format returns the format of the file (always "dcd") associated to the error
func (err Error) Format() string { return "dcd" }

//Critical returns true if the error is critical, false otherwise
func (err Error) Critical() bool { return err.critical }

const (
	TrajUnIni           = "Traj object uninitialized to read"
	ReadError           = "Error reading frame"
	UnableToOpen        = "Unable to open file"
	SecurityCheckFailed = "FailedSecurityCheck"
	WrongFormat         = "Wrong format in the DCD file or frame"
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

func (E lastFrameError) Format() string { return "dcd" }

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
