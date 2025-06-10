// +build !goyira
/*
 * dcd_write.go, part of gochem
 *
 * Copyright 2021 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
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
 *
 */

package dcd

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"runtime"

	v3 "github.com/rmera/gochem/v3"
)

//const mAXTITLE int32 = 80
//const rSCAL32BITS int32 = 1

//A writing buffer for DCD format
type WB []byte

//Writes w to the buffer
func (B WB) Write(w []byte) (int, error) {
	if len(B) < len(w) {
		return 0, fmt.Errorf("mismatched buffers")
	}
	for i := range w {
		B[i] = w[i]
	}
	return len(w), nil
}

//Container for an Charmm/NAMD binary trajectory file.
//opened for writing
type DCDWObj struct {
	natoms     int32
	buffSize   int
	readLast   bool //Have we read the last frame?
	writable   bool //Is it ready to be written on
	filename   string
	charmm     bool //Charmm traj?
	extrablock bool
	fourdim    bool
	frames     int32
	new        bool     //Still no frame written to it it?
	fixed      int32    //Fixed atoms (not supported)
	dcd        *os.File //The DCD file
	dcdFields  [][]float32
	concBuffer [][][]float32
	endian     binary.ByteOrder
}

//New writer initializes a DCD trajectory for writing.
func NewWriter(filename string, natoms int) (*DCDWObj, error) {
	traj := new(DCDWObj)
	traj.natoms = int32(natoms)
	if err := traj.initWrite(filename); err != nil {
		return nil, errDecorate(err, "New")
	}
	//	traj.dcdFields = make([][]float32, 3, 3)
	//	traj.dcdFields[0] = make([]float32, natoms, natoms)
	//	traj.dcdFields[1] = make([]float32, natoms, natoms)
	//	traj.dcdFields[2] = make([]float32, natoms, natoms)

	return traj, nil

}

func (D *DCDWObj) Close() {
	if !D.writable {
		return
	}
	D.dcd.Close()
	D.writable = false

}

//InitRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid.
//It support big and little endianness, charmm or (namd>=2.1) and no
//fixed atoms.
func (D *DCDWObj) initWrite(name string) error {
	wrapbinerr := func(err error) error {
		return Error{err.Error(), D.filename, []string{"binary.Write", "initWrite"}, true}
	}

	D.endian = binary.LittleEndian
	var err error
	D.dcd, err = os.Create(name)
	if err != nil {
		return Error{err.Error(), D.filename, []string{"os.Open", "initWrite"}, true}
	}
	if err := binary.Write(D.dcd, D.endian, int32(84)); err != nil {
		return wrapbinerr(err)
	}
	//For some reason, we have to write this magic number.
	magic := []byte("CORD")
	if err := binary.Write(D.dcd, D.endian, magic); err != nil {
		return wrapbinerr(err)
	}
	//The frames in the file go here. No frames written yet, but will update this part after every write.
	if err := binary.Write(D.dcd, D.endian, int32(0)); err != nil {
		return wrapbinerr(err)
	}
	//Initial time
	if err := binary.Write(D.dcd, D.endian, int32(0)); err != nil {
		return wrapbinerr(err)
	}
	//step interval (nsavc)
	if err := binary.Write(D.dcd, D.endian, int32(1)); err != nil {
		return wrapbinerr(err)
	}
	//5 zeros plus natom-rfreat
	for i := 0; i < 6; i++ {
		if err := binary.Write(D.dcd, D.endian, int32(0)); err != nil {
			return wrapbinerr(err)
		}

	}
	//delta time
	if err := binary.Write(D.dcd, D.endian, float32(1)); err != nil {
		return wrapbinerr(err)
	}
	//No unit cell
	if err := binary.Write(D.dcd, D.endian, int32(0)); err != nil {
		return wrapbinerr(err)
	}
	//8 zeros for charmm
	for i := 0; i < 8; i++ {
		if err := binary.Write(D.dcd, D.endian, int32(0)); err != nil {
			return wrapbinerr(err)
		}
	}
	//charmm version, let's say, 24
	if err := binary.Write(D.dcd, D.endian, int32(24)); err != nil {
		return wrapbinerr(err)
	}

	//don't ask me why
	if err := binary.Write(D.dcd, D.endian, int32(84)); err != nil {
		return wrapbinerr(err)
	}
	//same as above
	if err := binary.Write(D.dcd, D.endian, int32(244)); err != nil {
		return wrapbinerr(err)
	}

	//how many units of mAXTITLE does the title have?
	var ntitle int32 = 2 //just a dummy title.
	if err := binary.Write(D.dcd, D.endian, ntitle); err != nil {
		return wrapbinerr(err)
	}
	title := make([]byte, 2*mAXTITLE, 2*mAXTITLE)
	//Not a very good title, I suppose
	for j := range title {
		title[j] = byte('l')
	}
	title[len(title)-1] = byte('\000') //null-ended
	if err := binary.Write(D.dcd, D.endian, title); err != nil {
		return wrapbinerr(err)
	}
	//no idea
	if err := binary.Write(D.dcd, D.endian, int32(244)); err != nil {
		return wrapbinerr(err)
	}
	//no idea
	if err := binary.Write(D.dcd, D.endian, int32(4)); err != nil {
		return wrapbinerr(err)
	}
	//ok, this is important, the number of atoms in each snapshot
	//We should have got the number of atoms when we created the object. Will check just in case.
	//if it's zero it means it hasn't been set.
	if D.natoms == 0 {
		return Error{"Trajectory not initialized correctly, the number of atoms is set to zero!", D.filename, []string{"initWrite"}, true}
	}
	if err := binary.Write(D.dcd, D.endian, D.natoms); err != nil {
		return wrapbinerr(err)
	}
	//same as above
	if err := binary.Write(D.dcd, D.endian, int32(4)); err != nil {
		return wrapbinerr(err)
	}
	runtime.SetFinalizer(D, func(D *DCDWObj) {
		D.dcd.Close()
	})
	D.writable = true
	D.new = true //nothing read yet
	return nil   //nothing else to do
}

//WNext rites the next frame to the trajectory.
//the box isn't actually used, so far. It's only there for compatibility.
func (D *DCDWObj) WNext(towrite *v3.Matrix, box ...[]float64) error {
	if !D.writable {
		return Error{TrajUnIni, D.filename, []string{"WNext"}, true}
	}
	if towrite == nil {
		return Error{"got nil coordinates", D.filename, []string{"WNext"}, true}

	}
	if int32(towrite.NVecs()) != D.natoms {
		return Error{"Coordinates don't match the trajectory size", D.filename, []string{"WNext"}, true}
	}
	if D.dcdFields == nil {
		D.dcdFields = make([][]float32, 3, 3)
		D.dcdFields[0] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[1] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[2] = make([]float32, int(D.natoms), int(D.natoms))
	}
	//This is easier to write to the dcd
	for i := 0; i < int(D.natoms); i++ {
		k := i - i*(i/int(D.natoms)) //magic xD
		D.dcdFields[0][k] = float32(towrite.At(k, 0))
		D.dcdFields[1][k] = float32(towrite.At(k, 1))
		D.dcdFields[2][k] = float32(towrite.At(k, 2))
	}
	D.wnextRaw(D.dcdFields)
	D.frames++
	D.updateFrames()
	return nil
}

//Next Reads the next frame in a XtcObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DCDWObj) wnextRaw(blocks [][]float32) error {
	if len(blocks[0]) != int(D.natoms) || len(blocks[1]) != int(D.natoms) || len(blocks[2]) != int(D.natoms) {
		return Error{NotEnoughSpace, D.filename, []string{"nextRaw"}, true}
	}
	D.new = false
	var blocksize int32 = int32(len(blocks[0])) * 4 //the "4" is because the size is required in bytes, apparently.
	//now get the coords, each as a slice of float32
	//X
	if err := binary.Write(D.dcd, D.endian, blocksize); err != nil {
		return err
	}

	err := D.writeFloat32Block(blocks[0])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//Y
	//Collect the size first, then the rest
	if err := binary.Write(D.dcd, D.endian, blocksize); err != nil {
		return err
	}
	err = D.writeFloat32Block(blocks[1])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//Z
	if err := binary.Write(D.dcd, D.endian, blocksize); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Write", "wnextRaw"}, true}
	}
	err = D.writeFloat32Block(blocks[2])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//	fmt.Println("Z", blocks[2])
	//we skip the 4-D values if they exist. Apparently this is not present in the
	//last snapshot, so we use an EOF here to signal that we have read the last snapshot.
	return nil

}

//Writes a block of float32s to the file, adding its size
func (D *DCDWObj) writeFloat32Block(block []float32) error {
	var blocksize int32 = int32(len(block)) * 4
	if err := binary.Write(D.dcd, D.endian, block); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Write", "writeFloat32Block"}, true}

	}
	if err := binary.Write(D.dcd, D.endian, blocksize); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Write", "writeFloat32Block"}, true}
	}
	return nil
}

//DCD is silly enough to require the number of frames at the begining.
func (D *DCDWObj) updateFrames() error {
	currentoffset, err := D.dcd.Seek(0, io.SeekCurrent) //we'll need it to go back
	if err != nil {
		return Error{err.Error(), D.filename, []string{"dcd.Seek", "updateFrame"}, true}
	}
	//now we go to the begining of the file
	_, err = D.dcd.Seek(0, io.SeekStart)
	if err != nil {
		return Error{err.Error(), D.filename, []string{"dcd.Seek", "updateFrame"}, true}
	}

	wrapbinerr := func(err error) error {
		return Error{err.Error(), D.filename, []string{"binary.Write", "updateFrame"}, true}
	}
	//I could try to get directly to the part we need to write, but it's so close to the begining that I think it's best to just write a couple
	//of unnecesary numbers.
	if err := binary.Write(D.dcd, D.endian, int32(84)); err != nil {
		return wrapbinerr(err)
	}
	//For some reason, we have to write this magic number.
	magic := []byte("CORD")
	if err := binary.Write(D.dcd, D.endian, magic); err != nil {
		return wrapbinerr(err)
	}
	if err := binary.Write(D.dcd, D.endian, int32(D.frames)); err != nil {
		return wrapbinerr(err)
	}
	//we go back to the part of the file we were reading
	_, err = D.dcd.Seek(currentoffset, io.SeekStart)
	if err != nil {
		return Error{err.Error(), D.filename, []string{"dcd.Seek", "updateFrame"}, true}
	}

	return nil

}
