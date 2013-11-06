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

import "fmt"
import "encoding/binary"
import "os"
import "bytes"
import "runtime"
import "github.com/rmera/gochem"

const MAXTITLE int32 = 80
const RSCAL32BITS int32 = 1

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

func NewDCD(filename string) (*DCDObj, error) {
	traj := new(DCDObj)
	if err := traj.initRead(filename); err != nil {
		return nil, err
	}
	traj.dcdFields = make([][]float32, 3, 3)
	traj.dcdFields[0] = make([]float32, int(traj.natoms), int(traj.natoms))
	traj.dcdFields[1] = make([]float32, int(traj.natoms), int(traj.natoms))
	traj.dcdFields[2] = make([]float32, int(traj.natoms), int(traj.natoms))
	traj.concBuffer = append(traj.concBuffer, traj.dcdFields)
	return traj, nil

}

//Returns true if the object is ready to be read from
//false otherwise. It doesnt guarantee that there is something
//to read.
//true or false depending on whether D is ready to read
//snapshots from it.
func (D *DCDObj) Readable() bool {
	return D.readable
}

//InitRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid.
//It support big and little endianness, charmm or (namd>=2.1) and no
//fixed atoms.
func (D *DCDObj) initRead(name string) error {
	rec_scale := RSCAL32BITS //At least for now we will not support anything else.
	D.endian = binary.LittleEndian
	_ = rec_scale
	NB := bytes.NewBuffer //shortness sake
	var err error
	D.dcd, err = os.Open(name)
	if err != nil {
		return err
	}
	var check int32
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return err
	}
	//For some reason the first thing we should read is an 84.
	//If this fails it means that the file is big endian.
	if check != 84 {
		D.endian = binary.BigEndian //
	}
	//Then the magic number "CORD", also for some unknown reason.
	magic := make([]byte, 4, 4)
	if err := binary.Read(D.dcd, D.endian, magic); err != nil {
		return err
	}
	if string(magic) != "CORD" {
		return fmt.Errorf("Wrong magic number")
	}

	//We first read a big chuck for random access.
	buf := make([]byte, 80, 80)
	if err := binary.Read(D.dcd, D.endian, buf); err != nil {
		return err
	}
	//X-plor sets this last int to zero, charmm sets it to its version number.
	//if we have a charmm file we get some additional flags.
	if err := binary.Read(NB(buf[76:]), D.endian, &check); err != nil {
		return err
	}
	if check != 0 {
		//		fmt.Println("CHARMM!!!") //////77
		D.charmm = true
		if err := binary.Read(NB(buf[40:]), D.endian, &check); err != nil {
			return err
		}
		if check != 0 {
			//			fmt.Println("block", check) ///////////
			D.extrablock = true
		}
		if err := binary.Read(NB(buf[40:]), D.endian, &check); err != nil {
			return err
		}
		if check == 1 {
			//			fmt.Println("4-dim", check) ///////////
			D.fourdim = true
		}

	} else {
		return fmt.Errorf("X-plor DCD not supported")
	}
	if err := binary.Read(NB(buf[32:]), D.endian, &D.fixed); err != nil {
		return err
	}
	//	fmt.Println("fixed", D.fixed)
	var delta float32 //This should work only on Charmm and namd >=2.1
	if err := binary.Read(NB(buf[36:]), D.endian, &delta); err != nil {
		return err
	}
	//	fmt.Println("delta:", delta)///////////////////////////////////////

	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return err
	}
	if check != 84 {
		return fmt.Errorf("Wrong DCD format")
	}
	var input_int int32
	if err := binary.Read(D.dcd, D.endian, &input_int); err != nil {
		return err
	}
	//how many units of MAXTITLE does the title have?
	var ntitle int32
	if err := binary.Read(D.dcd, D.endian, &ntitle); err != nil {
		return err
	}
	title := make([]byte, MAXTITLE*ntitle, MAXTITLE*ntitle)
	if err := binary.Read(D.dcd, D.endian, title); err != nil {
		return err
	}
	//	fmt.Println("Title:", string(title))///////////////////////////////////////
	if err := binary.Read(D.dcd, D.endian, &input_int); err != nil {
		return err
	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return err
	}
	if check != 4 { //one must read a 4 before the natoms
		return fmt.Errorf("Wrong format in DCD")
	}
	if err := binary.Read(D.dcd, D.endian, &D.natoms); err != nil {
		return err
	}
	//	fmt.Println("natoms", D.natoms)
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return err
	}
	if check != 4 { //and one more 4
		return fmt.Errorf("DCD has wrong format")
	}
	if D.fixed == 0 {
		runtime.SetFinalizer(D, func(D *DCDObj) {
			D.dcd.Close()
		})
		D.readable = true
		return nil //nothing else to do
	}
	D.new = true //nothing read yet
	return fmt.Errorf("Fixed atoms not supported")

}

//Next Reads the next frame in a DcDObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DCDObj) Next(keep *chem.VecMatrix) error {
	if !D.readable {
		return fmt.Errorf("Not readable")
	}
	if D.dcdFields == nil {
		D.dcdFields = make([][]float32, 3, 3)
		D.dcdFields[0] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[1] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[2] = make([]float32, int(D.natoms), int(D.natoms))
	}
	if err := D.nextRaw(D.dcdFields); err != nil {
		return D.eOF2NoMoreFrames(err)
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
		return fmt.Errorf("Not enough space in passed blocks")
	}
	D.new = false
	if D.readLast {
		D.readable = false
		return fmt.Errorf("No more frames")
	}
	//if there is an extra block we just skip it.
	//Sadly, even when there is an extra block, it is not present in all
	//snapshots for some trajectories, so we must use the block size to see if
	//there is an extra block or if the X block starts inmediately
	var blocksize int32
	if D.extrablock {
		if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
			return err
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
			return err
		}
	}
	err := D.readFloat32Block(blocksize, blocks[0])
	if err != nil {
		return err
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
		return err
	}
	//	fmt.Println("Y", blocks[1])
	//Z
	if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
		return err
	}
	err = D.readFloat32Block(blocksize, blocks[2])
	if err != nil {
		return err
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
				return err
			}
		}
		if !D.readLast {
			if _, err := D.readByteBlock(blocksize); err != nil {
				return err
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
		return err
	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return err
	}
	if check != blocksize {
		return fmt.Errorf("Wrong format in DCD snapshot")
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
		return nil, err
	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return nil, err
	}
	if check != blocksize {
		return nil, fmt.Errorf("Failed security check")
	}
	return block, nil
}

//Natoms returns the number of atoms per frame in the XtcObj.
//XtcObj must be initialized. 0 means an uninitialized object.
func (D *DCDObj) Len() int {
	return int(D.natoms)
}

func (D *DCDObj) eOF2NoMoreFrames(err error) error {
	if err == nil {
		return nil
	}
	if err.Error() == "EOF" {
		return fmt.Errorf("No more frames")
	}
	return err
}

func (D *DCDObj) setConcBuffer(batchsize int) error {
	l := D.buffSize
	if l == batchsize {
		return nil
	} else if l > batchsize {
		for i := batchsize; i < l; i++ {
			for j, _ := range D.concBuffer[i] {
				D.concBuffer[i][j] = nil
			}
			D.concBuffer[i] = nil //no idea if this actually works
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

/*NextConc takes a slice of bools and reads as many frames as elements the list has
form the trajectory. The frames are discarted if the corresponding elemetn of the slice
* is false. The function returns a slice of channels through each of each of which
* a *matrix.DenseMatrix will be transmited*/
func (D *DCDObj) NextConc(frames []*chem.VecMatrix) ([]chan *chem.VecMatrix, error) {
	if !D.Readable() {
		return nil, fmt.Errorf("Traj object uninitialized to read")
	}
	framechans := make([]chan *chem.VecMatrix, len(frames)) //the slice of chans that will be returned
	if D.buffSize < len(frames) {
		D.setConcBuffer(len(frames))
	}
	for key, _ := range frames {
		DFields := D.concBuffer[key]
		if err := D.nextRaw(DFields); err != nil {
			return nil, D.eOF2NoMoreFrames(err)
		}
		//We have to test for used twice to allow allocating for goCoords
		//When the buffer is not going to be used.
		if frames[key] == nil {
			framechans[key] = nil //ignored frame
			continue
		}
		framechans[key] = make(chan *chem.VecMatrix)
		//Now the parallel part
		go func(natoms int, DFields [][]float32, keep *chem.VecMatrix, pipe chan *chem.VecMatrix) {
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
