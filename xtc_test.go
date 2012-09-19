// +build gromacs 

/*
 * untitled.go
 * 
 * Copyright 2012 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by the Free Software Foundation; either version 2.1 of the 
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General 
 * Public License along with this program.  If not, see 
 * <http://www.gnu.org/licenses/>.
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
import "testing"
import "github.com/skelterjohn/go.matrix"

/*TestXtc reads the frames of the test xtc file using the
 * "interactive" or "low level" functions, i.e. one frame at a time
 * It prints the firs 2 coordinates of each frame and the number of 
 * read frames at the end.*/
func TestXtc(Te *testing.T) {
	fmt.Println("Fist test!")
	name:="test/test.xtc"
	traj:=new(XtcObj)
	err:=traj.InitRead(name)
	if err!=nil{
		Te.Error(err)
		}
	i:=0
	for ;;i++{
		coords,err:=traj.Next(true)
		if err!=nil && err.Error()!="No more frames"{
			Te.Error(err)
			break
			}else if err==nil{
			fmt.Println(coords.GetRowVector(2))
			}else{
			break	
			}
		}
	fmt.Println("Over! frames read:", i)
	}

/*TestFrameXtc reads the frames of the test xtc file from the first to
 * the forth frame skipping one frame for each read one. It uses the
 * "high level" function. It prints the frames read twince, and the
 * coordinates of the forth atom of the last read frame,*/	
func TestFrameXtc(Te *testing.T) {
	fmt.Println("Second test!")
	name:="test/test.xtc"
	traj:=new(XtcObj)
	err:=traj.InitRead(name)
	if err!=nil{
		Te.Error(err)
		}
	Coords,read,err:=traj.ManyFrames(0, 5, 1)
	if err!=nil{
		Te.Error(err)
		}
	fmt.Println(len(Coords),read,Coords[read-1].GetRowVector(4))
	}

func TestFrameXtcConc(Te *testing.T){
	fmt.Println("Third test!")
	name:="test/test.xtc"
	traj:=new(XtcObj)
	err:=traj.InitRead(name)
	if err!=nil{
		Te.Error(err)
		}
	i:=0
	for ;;i++{
		_=matrix.Zeros(3,3)
		frames:=[]bool{true,true}
		coordschans,err:=traj.NextConc(frames)
		if err!=nil && err.Error()!="No more frames"{
			Te.Error(err)
			break
			}else if err==nil{
			for _,pipe:=range(coordschans){
				fmt.Println((<-pipe).GetRowVector(2))
				}
			
			}else{
			break	
			}
		}
	}	
		/*
			results:=make([]chan *matrix.DenseMatrix,len(coordchans))
			for key,channel:=range(coordchans){
				results=append(results,make(chan *matrix.DenseMatrix,1))
				go func(channelin,channelout chan *matrix.DenseMatrix,current int){
					if channelin!=nil{
						matrix:=<-channelin
						vector:=matrix.GetRowVector(2)
					//	fmt.Println(current, vector)
						channelout<-vector
						}else{
			//			fmt.Println("This frame was dropped",current)	
						channelout<-nil
						}
					}(channel,results[key],i+(2*i)+key)
				}
			for _,j:=range(results){
				fmt.Println(<-j)
				}
			}else{
			break	
			}
		}
	fmt.Println("Over! frames read:", i)
	}

/*						channelout<-vector
						}else{
						channelout<-nil	
						}
					}(channel,results[key])
				}
			for _,j:=range(results){
				fmt.Println(<-j)
				}
*/
