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
	"os"
	"runtime"

	v3 "github.com/rmera/gochem/v3"
)

//const mAXTITLE int32 = 80
//const rSCAL32BITS int32 = 1

//Container for an Charmm/NAMD binary trajectory file.
type DCDWObj struct {
	natoms     int32
	buffSize   int
	readLast   bool //Have we read the last frame?
	writable   bool //Is it ready to be written on
	filename   string
	charmm     bool //Charmm traj?
	extrablock bool
	fourdim    bool
	new        bool     //Still no frame written to it it?
	fixed      int32    //Fixed atoms (not supported)
	dcd        *os.File //The DCD file
	dcdFields  [][]float32
	concBuffer [][][]float32
	endian     binary.ByteOrder
}

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

//InitRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid.
//It support big and little endianness, charmm or (namd>=2.1) and no
//fixed atoms.
func (D *DCDWObj) initWrite(name string) error {
	wrapbinerr := func(err error) error {
		return Error{err.Error(), D.filename, []string{"binary.Write", "initWrite"}, true}
	}

	D.endian = binary.LittleEndian
	NB := bytes.NewBuffer //shortness sake
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
	//We first read a big chuck for random access.
	buf := make([]byte, 80, 80)
	//X-plor sets this last int to zero, charmm sets it to its version number.
	//if we have a charmm file we get some additional flags.
	if err := binary.Write(NB(buf[76:]), D.endian, int32(2)); err != nil {
		return wrapbinerr(err)
	}
	D.charmm = true
	//no extra block of fourth dimmension
	if err := binary.Write(NB(buf[40:]), D.endian, int32(0)); err != nil {
		return wrapbinerr(err)
	}
	//here go the fixed atoms, I just set it to 0
	if err := binary.Write(NB(buf[32:]), D.endian, D.fixed); err != nil {
		return Error{err.Error(), D.filename, []string{"initWrite"}, true}

	}
	var delta float32 = 1.0 //This should work only on Charmm and namd >=2.1
	//no idea what delta is. Right now, I just write 0.
	//I really have to check this with the format!
	if err := binary.Write(NB(buf[36:]), D.endian, delta); err != nil {
		return wrapbinerr(err)
	}
	//we now write the buffer we had prepared to the file.
	binary.Write(D.dcd, D.endian, buf)
	//Again, no idea why the number 84 goes here. But it does.
	if err := binary.Write(D.dcd, D.endian, int32(84)); err != nil {
		return Error{err.Error(), D.filename, []string{"initWrite"}, true}
	}
	//another number that I seem to have ignored when reading. Will just write 0
	if err := binary.Write(D.dcd, D.endian, int32(0)); err != nil {
		return wrapbinerr(err)
	}
	//how many units of mAXTITLE does the title have?
	var ntitle int32 = 1 //just a dummy title of the smallest size possible.
	if err := binary.Write(D.dcd, D.endian, ntitle); err != nil {
		return wrapbinerr(err)
	}
	title := make([]byte, mAXTITLE, mAXTITLE)
	if err := binary.Write(D.dcd, D.endian, title); err != nil {
		return wrapbinerr(err)
	}
	//another unknown value, agains, I'm just writing 0 here, because I don't know what this is.
	if err := binary.Write(D.dcd, D.endian, int32(0)); err != nil {
		return wrapbinerr(err)

	}
	//For some reason, there must be a 4 before the natoms.
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
	//one more fore for some reason.
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

//Writes the next frame to the trajectory.
func (D *DCDWObj) WNext(towrite *v3.Matrix) error {
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
		k := i - i*(i/int(D.natoms))
		D.dcdFields[0][k] = float32(towrite.At(k, 0))
		D.dcdFields[1][k] = float32(towrite.At(k, 1))
		D.dcdFields[2][k] = float32(towrite.At(k, 2))
	}
	D.wnextRaw(D.dcdFields)
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
	var blocksize int32 = int32(len(blocks[0]))
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
	var blocksize int32 = int32(len(block))
	if err := binary.Write(D.dcd, D.endian, block); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Write", "writeFloat32Block"}, true}

	}
	if err := binary.Write(D.dcd, D.endian, blocksize); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Write", "writeFloat32Block"}, true}
	}
	return nil
}
