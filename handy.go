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
	//	fmt.Println(reflect.TypeOf(set).Kind().String())
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
	
	
	
	
	
func RotateAbout(angle float64,axis, coordsorig *matrix.DenseMatrix) (*matrix.DenseMatrix,error){
	coords:=coordsorig.Copy()
	translation:=coords.GetRowVector(0).Copy()
	_=translation
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
