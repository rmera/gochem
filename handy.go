/*
 * files.go, part of gochem.
 * 
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 * 
 * 
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 * 
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

	
package chem



import "reflect"
import "github.com/skelterjohn/go.matrix"
import "fmt"



func Deg2Rad(f float64)float64{
	return f*0.0174533
	}

func Rad2Deg(f float64)float64{
	return f/0.0174533
	}


/*IsIn returns the position of test in the slice set, or
 * -1 if test is not present in set. Panics if set is not a slice*/
func IsIn(test interface{}, set interface{}) (int){
	vset:=reflect.ValueOf(set)
	if reflect.TypeOf(set).Kind().String()!="slice"{
		panic("IsIn function needs a slice as second argument!")
		}
	if vset.Len()<0{
		return 1
		}
	for i:=0;i<vset.Len();i++{
		vcomp:=vset.Index(i)
		comp:=vcomp.Interface()
		if reflect.DeepEqual(test, comp){
			return i
			}
		}
	return -1
	}
	
	
//Ones mass returns a column matrix with lenght rosw. 
//This matrix can be used as a dummy mass matrix
//for geometric calculations.
func OnesMass(lenght int)*matrix.DenseMatrix{
	return matrix.Ones(lenght,1)
	}


	
//Super determines the best rotation and translations to superimpose the coords in test
//listed in testlst on te atoms of molecule templa, frame frametempla, listed in templalst. 
//It applies those rotation and translations to the whole frame frametest of molecule test, in palce. 
//testlst and templalst must have the same number of elements.
func Super(test, templa *matrix.DenseMatrix, testlst, templalst []int) (*matrix.DenseMatrix, error){
	ctest:=test
	if len(testlst)!=0{
		ctest=SomeRows(test,testlst)
		}
	ctempla:=templa
	if len(templalst)!=0{
		ctempla=SomeRows(templa,templalst,)
		}
	if ctempla.Rows()!=ctest.Rows(){
		return nil, fmt.Errorf("Mismatched template and test atom numbers: %d, %d", ctempla.Rows(),ctest.Rows())
		}
	_,rotation,trans1,trans2,err1:=GetSuper(ctest,ctempla)
	if err1!=nil{
		return nil,err1
		}
	err1=AddRow(test,trans1)
	test=matrix.ParallelProduct(test,rotation)
	err2:=AddRow(test,trans2)
	if err1 != nil || err2!=nil{
		return nil,fmt.Errorf("Unexpected error when aplying superposition")
		}
	return test,nil
	}



//Rotate about rotates the coordinates in coordsorig around by angle radians around the axis 
//given by the vector axis. It returns the rotated coordsorig, since the original is not affected.
func RotateAbout(coordsorig, ax1, ax2 *matrix.DenseMatrix,angle float64) (*matrix.DenseMatrix,error){
	coords:=coordsorig.Copy()
	translation:=ax1.Copy()
	axis:=ax1.Copy()
	axis.Subtract(ax2)//now it became the rotation axis
	err:=SubRow(coords,translation)
	if err!=nil{
		return nil, err
		}
	Zswitch:=GetSwitchZ(axis)
	coords=matrix.ParallelProduct(coords,Zswitch) //rotated
	Zrot,err:=GetRotateAroundZ(angle)
	if err!=nil{
		return nil, err
		}
	RevZ,err:=Zswitch.Inverse()
	if err!=nil{
		return nil, err
		}
	coords=matrix.ParallelProduct(coords,Zrot) //rotated
	coords=matrix.ParallelProduct(coords,RevZ)
	err=AddRow(coords,translation)
	if err!=nil{
		return nil, err
		}
	return coords, err
	}


//Corrupted is a convenience function to check that a reference and a trajectory have the same number of atoms
func Corrupted(R Ref,X Traj)error{
	if X.Len()!=R.Len(){
		return fmt.Errorf("Mismatched number of atoms/coordinates")
		}
	return nil
	}

