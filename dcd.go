////////// +build  D.dcd

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
	readLast bool  //Have we read the last frame?
	readable bool  //Is it ready to be read
	filename string
	charmm bool     //Charmm traj?
	extrablock bool
	fourdim bool
	new bool     //Still no frame read from it?
	fixed int32   //Fixed atoms (not supported)
	dcd *os.File //The DCD file
	dcdFields [][]float32
	
	
}

//Returns true if the object is ready to be read from
//false otherwise. It doesnt guarantee that there is something
//to read.
//true or false depending on whether D is ready to read
//snapshots from it.
func (D *DcdObj) Readable() bool{
		return D.readable
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
		D.readable=true
		return nil //nothing else to do
		}
	D.new=true //nothing read yet
	return fmt.Errorf("Fixed atoms not supported")


}

//Next Reads the next frame in a DcDObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DcdObj) Next(keep bool) (*matrix.DenseMatrix, error) {
	if !D.readable{
		return nil, fmt.Errorf("Not readable")
		}
	if D.dcdFields==nil{
		D.dcdFields=make([][]float32,3,3)
		D.dcdFields[0]=make([]float32,int(D.natoms),int(D.natoms))
		D.dcdFields[1]=make([]float32,int(D.natoms),int(D.natoms))
		D.dcdFields[2]=make([]float32,int(D.natoms),int(D.natoms))
		}
	if err:=D.nextRaw(D.dcdFields);err!=nil{
		return nil, err
		}
	if !keep{
		return nil, nil
		}
	outBlock:=make([]float64,int(D.natoms)*3,int(D.natoms)*3)
	for i:=0;i<int(D.natoms);i++{
		j:=i+(i*2)
		k:=i-i*(i/int(D.natoms))
		outBlock[j]=float64(D.dcdFields[0][k])
		outBlock[j+1]=float64(D.dcdFields[1][k])
		outBlock[j+2]=float64(D.dcdFields[2][k])
	}
	final:=matrix.MakeDenseMatrix(outBlock,int(D.natoms),3)
	fmt.Print(final)/////////7
	return final, nil
	
}



//Next Reads the next frame in a XtcObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DcdObj) nextRaw(blocks [][]float32) (error) {
	if len(blocks[0])!=int(D.natoms) || len(blocks[1])!=int(D.natoms) || len(blocks[2])!=int(D.natoms){
		return fmt.Errorf("Not enough space in passed blocks")	
	}
	D.new=false
	if D.readLast{
		D.readable=false
		return fmt.Errorf("No more frames")
		}
	
	//if there is an extra block we just skip it.
	//Sadly, even when there is an extra block, it is not present in all
	//snapshots for some trajectories, so we must use the block size to see if
	//there is an extra block or if the X block starts inmediately
	var blocksize int32
	if D.extrablock{
		if err:=binary.Read(D.dcd,binary.LittleEndian,&blocksize);err!=nil{
			return err
		}
		//If the blocksize is 4*natoms it means that the block is not an 
		//extra block, but the X coordinates, and thus we must skip the following
		if blocksize!=D.natoms*4{ 
			if _,err:=D.readByteBlock(blocksize);err!=nil{ 
				return err
			}
			blocksize=0
		}
	}
	//now get the coords, each as a slice of float32
	//X
	//we collect the X block size again only if it has not been collected before
	if blocksize==0{
			if err:=binary.Read(D.dcd,binary.LittleEndian,&blocksize);err!=nil{
				return err
			}
		}
    err:=D.readFloat32Block(blocksize,blocks[0])
	if err!=nil{
		if err.Error()=="EOF"{
			return fmt.Errorf("No more frames")
			}
		return err
		}
//	fmt.Println("X", len(xblock)) //, xblock)
	fmt.Println("X", blocks[0])  ///////////////////////////////
	//Y
	//Collect the size first, then the rest
	if err:=binary.Read(D.dcd,binary.LittleEndian,&blocksize);err!=nil{
		return  err
		}
	err=D.readFloat32Block(blocksize,blocks[1])
	if err!=nil{
		return err
		}
	fmt.Println("Y", blocks[1]) 
	//Z
	if err:=binary.Read(D.dcd,binary.LittleEndian,&blocksize);err!=nil{
		return err
		}
	err=D.readFloat32Block(blocksize,blocks[2])
	if err!=nil{
		return err
		}
	fmt.Println("Z", blocks[2]) 
	//we skip the 4-D values if they exist. Apparently this is not present in the 
	//last snapshot, so we use an EOF here to signal that we have read the last snapshot.
	if D.charmm && D.fourdim{
		if err:=binary.Read(D.dcd,binary.LittleEndian,&blocksize);err!=nil{
			if err.Error()=="EOF"{
				D.readLast=true
			//	fmt.Println("LAST!")
			}else{
				return err
			}
		}
		if !D.readLast{
			if _,err:=D.readByteBlock(blocksize);err!=nil{
				return err
			}
		}
	}
	return  nil
	
}

//Queries the size of a block, and reads its contents into block, which must have the 
//appropiate size.
func (D *DcdObj)readFloat32Block(blocksize int32, block []float32)(error) {
	var check int32
	fmt.Println("blockf",blocksize)	
	if err:=binary.Read(D.dcd,binary.LittleEndian,block);err!=nil{
		return err
	}
	if err:=binary.Read(D.dcd,binary.LittleEndian,&check);err!=nil{
		return err
	} 	
	if check!=blocksize{
		return fmt.Errorf("Wrong format in DCD snapshot")
	}
	return nil	
}


//Queries the size of a block, make a slice of a quarter of that size
//and reads that ammount of float32. This function is used for the
//
func (D *DcdObj)readByteBlock(blocksize int32)([]byte,error) {
	var check int32
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
