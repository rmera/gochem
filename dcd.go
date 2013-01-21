////////// +build  D.dcd

/*
 * D.dcd.go, part of gochem
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


import "fmt"
import "github.com/skelterjohn/go.matrix"
import "encoding/binary"
import "os"
import "bytes"
import "runtime"

const MAXTITLE int32 = 80
const RSCAL32BITS int32 = 1

//Container for an Charmm/NAMD binary trajectory file.
type DcdObj struct {
	natoms   int32
	filename string
	charmm bool
	extrablock bool
	fourdim bool
	fixed int32
	dcd *os.File
	
}

//Returns true if the object is ready to be read from
//false otherwise. It doesnt guarantee that there is something
//to read.
func (D *DcdObj) Readable() bool {
	return true
	
}

//InitRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid.
//It only support little-endianness, charmm/namd>=2.1 and no
//fixed atoms.
func (D *DcdObj) InitRead(name string) error {
	rec_scale:=RSCAL32BITS //At least for now we will not support anything else.
	_=rec_scale
	NB:=bytes.NewBuffer //shortness sake
	var err error
	D.dcd, err=os.Open(name)
	if err!=nil{
		return err
	}
	var check int32
	if err:=binary.Read(D.dcd,binary.LittleEndian,&check);err!=nil{
		return err
	} 
	//For some reason the first thing we should read is an 84.
	if check!=84{
		return fmt.Errorf("Endianness probably wrong")
	}
	//Then the magic number "CORD", also for some unknown reason.
	magic:=make([]byte,4,4)
	if err:=binary.Read(D.dcd,binary.LittleEndian,magic);err!=nil{
		return err
	} 
	if string(magic)!="CORD"{
		return fmt.Errorf("Wrong magic number")
	}
	
	//We first read a big chuck for random access.
	buf:=make([]byte,80,80)
	if err:=binary.Read(D.dcd,binary.LittleEndian,buf);err!=nil{
		return err
	} 	
	//X-plor sets this last int to zero, charmm sets it to its version number.
	//if we have a charmm file we get some additional flags.
	if err:=binary.Read(NB(buf[76:]),binary.LittleEndian,&check);err!=nil{
		return err
	} 	
	if check!=0{
		fmt.Println("CHARMM!!!") //////77
		D.charmm=true
		if err:=binary.Read(NB(buf[40:]),binary.LittleEndian,&check);err!=nil{
			return err
		}
		if check!=0{
			fmt.Println("block", check) ///////////
			D.extrablock=true
		}
		if err:=binary.Read(NB(buf[40:]),binary.LittleEndian,&check);err!=nil{
			return err
		}
		if check==1{
			fmt.Println("4-dim", check) ///////////
			D.fourdim=true
		}		
		
	}else{
		return fmt.Errorf("X-plor DCD not supported")
	}
	if err:=binary.Read(NB(buf[32:]),binary.LittleEndian,&D.fixed);err!=nil{
		return err
	}
	fmt.Println("fixed", D.fixed) 		
	var delta float32 //This should work only on Charmm and namd >=2.1
	if err:=binary.Read(NB(buf[36:]),binary.LittleEndian,&delta);err!=nil{
		return err
	}
	fmt.Println("delta:", delta)/////////////////////////////////////// 	

	if err:=binary.Read(D.dcd,binary.LittleEndian,&check);err!=nil{
		return err
	} 
	if check!=84{
		return fmt.Errorf("Wrong DCD format")
	}	
	var input_int int32
	if err:=binary.Read(D.dcd,binary.LittleEndian,&input_int);err!=nil{
		return err
	} 
	//how many units of MAXTITLE does the title have?
	var ntitle int32
	if err:=binary.Read(D.dcd,binary.LittleEndian,&ntitle);err!=nil{
		return err
	} 
	title:=make([]byte,MAXTITLE*ntitle,MAXTITLE*ntitle)
	if err:=binary.Read(D.dcd,binary.LittleEndian,title);err!=nil{
		return err
	}
	fmt.Println("Title:", string(title))/////////////////////////////////////// 
	if err:=binary.Read(D.dcd,binary.LittleEndian,&input_int);err!=nil{
		return err
	} 
	if err:=binary.Read(D.dcd,binary.LittleEndian,&check);err!=nil{
		return err
	} 
	if check!=4{  //one must read a 4 before the natoms
		return fmt.Errorf("Wrong format in DCD")
	}	
	if err:=binary.Read(D.dcd,binary.LittleEndian,&D.natoms);err!=nil{
		return err
	} 
	fmt.Println("natoms", D.natoms)
	if err:=binary.Read(D.dcd,binary.LittleEndian,&check);err!=nil{
		return err
	} 
	if check!=4{  //and one more 4
		return fmt.Errorf("Endianness probably wrong")
	}
	if D.fixed==0{
		runtime.SetFinalizer(D, func(D *DcdObj) {
		D.dcd.Close()
		})
		return nil //nothing else to do
		}
	return fmt.Errorf("Fixed atoms not supported")
	
/*	freeindexes:=make([]int32,D.natoms-D.fixed)
//	fixedcoords:=make([]float32,D.natoms*4-D.fixed)
	if err:=binary.Read(D.dcd,binary.LittleEndian,&input_int);err!=nil{
		return err
	}
//	fmt.Println(input_int)
	if input_int!=(D.natoms-D.fixed*4){
		return fmt.Errorf("Wrong format in DCD")
		} 	
	if err:=binary.Read(D.dcd,binary.LittleEndian,freeindexes);err!=nil{
		return err
	} 
	if err:=binary.Read(D.dcd,binary.LittleEndian,&input_int);err!=nil{
		return err
	}
//	fmt.Println(input_int)
	if input_int!=(D.natoms-D.fixed*4){
		return fmt.Errorf("Wrong format in DCD")
		} 

	
	return nil
*/
}

//Next Reads the next frame in a XtcObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DcdObj) Next(keep bool) (*matrix.DenseMatrix, error) {
	//if there is an extra block we just skip it.
	if D.charmm && D.extrablock && keep{
		if _,err:=D.readByteBlock();err!=nil{ 
			return nil, err
			}
	}
	//now get the coords, each as a slice of float32
	//X
	xblock, err:=D.readFloat32Block()
	if err!=nil{
		if err.Error()=="EOF"{
			return nil,fmt.Errorf("No more frames")
			}
		return nil,err
		}
	fmt.Println("X", len(xblock)) //, xblock)
	//Y
	yblock, err:=D.readFloat32Block()
	if err!=nil{
		return nil,err
		}
	fmt.Println("Y", len(yblock)) //, yblock)
	//Z
	zblock, err:=D.readFloat32Block()
	if err!=nil{
		return nil,err
		}
	fmt.Println("Z", len(zblock))//, zblock)
	//we skip the 4-D values if they exist
	if D.charmm && D.fourdim{
		if _,err:=D.readByteBlock();err!=nil{
			if err.Error()=="EOF"{
				return nil,fmt.Errorf("No more frames")
				} 
			return nil, err
			}
	}
	xlen:=len(xblock)
	ylen:=len(yblock)
	zlen:=len(zblock)
	if xlen!=ylen || ylen!=zlen || xlen!=int(D.natoms){
		return nil,fmt.Errorf("Wrong format in DCD snapshot") 
		}
	return nil, nil //Just drop the frame
	
}

//Queries the size of a block, make a slice of a quarter of that size
//and reads that ammount of float32. This function is used for the
//
func (D *DcdObj)readFloat32Block()([]float32,error) {
	var blocksize int32
	var check int32
	if err:=binary.Read(D.dcd,binary.LittleEndian,&blocksize);err!=nil{
		return nil, err
	} 
	fmt.Println("blockf",blocksize)	
	block:=make([]float32,blocksize/4,blocksize/4)
	if err:=binary.Read(D.dcd,binary.LittleEndian,block);err!=nil{
		return nil,err
	}
	if err:=binary.Read(D.dcd,binary.LittleEndian,&check);err!=nil{
		return nil,err
	} 	
	if check!=blocksize{
		return nil,fmt.Errorf("Wrong format in DCD snapshot")
	}
	return block, nil	
}


//Queries the size of a block, make a slice of a quarter of that size
//and reads that ammount of float32. This function is used for the
//
func (D *DcdObj)readByteBlock()([]byte,error) {
	var blocksize int32
	var check int32
	if err:=binary.Read(D.dcd,binary.LittleEndian,&blocksize);err!=nil{
		return nil, err
	} 
	fmt.Println("blockb",blocksize)	
	block:=make([]byte,blocksize,blocksize)
	if err:=binary.Read(D.dcd,binary.LittleEndian,block);err!=nil{
		return nil,err
	}
	if err:=binary.Read(D.dcd,binary.LittleEndian,&check);err!=nil{
		return nil,err
	} 	
	if check!=blocksize{
		return nil,fmt.Errorf("Failed security check")
	}
	return block, nil	
}





/*NextConc takes a slice of bools and reads as many frames as elements the list has
form the trajectory. The frames are discarted if the corresponding elemetn of the slice
* is false. The function returns a slice of channels through each of each of which 
* a *matrix.DenseMatrix will be transmited*/
func (D *DcdObj) NextConc(frames []bool) ([]chan *matrix.DenseMatrix, error) {
	return nil, nil
}

//Read frames from Traj from ini to end skipping skip frames between read. Returns a slice with coords of each frame
//the number of frames read and error or nil.
func (D *DcdObj) ManyFrames(ini, end, skip int) ([]*matrix.DenseMatrix, int, error) {
	return nil, 0, nil
}

//Natoms returns the number of atoms per frame in the XtcObj.
//XtcObj must be initialized. 0 means an uninitialized object.
func (D *DcdObj) Len() int {
	return 1
}


