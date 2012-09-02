// +build !gromacs

/*
 * general_test.go
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
 */

/*This provides some tests for the library, in the form of little functions 
 * that have practical applications*/

package chem


import "github.com/skelterjohn/go.matrix"
import "fmt"
import "testing"
import "strings"
import "strconv"




//TestChangeAxis reads the PDB 2c9v.pdb from the test directory, collects 
//The CA and CB of residue D124, and rotates the whole molecule such as the vector
//defined by these 2 atoms is aligned with the Z axis. The new molecule is written
//as 2c9v_aligned.pdb to the test folder.
func TestChangeAxis(Te *testing.T){
	var mol Molecule
	ats,coords,bfac,err:=PdbRead("test/2c9v.pdb",true)
	if err!=nil{
		Te.Error(err)
		}
	mol.Atoms=ats
	mol.Coords=coords
	mol.Bfactors=bfac
	orient_atoms:=[2]int{0,0}
	for index, atom:= range(mol.Atoms){
		if atom.Chain=='A' && atom.Molid==124{
			if atom.Name=="CA"{
				orient_atoms[0]=index
				}else if atom.Name=="CB"{
				orient_atoms[1]=index	
				}
			}
		}
	ov1:=mol.Coord(orient_atoms[0], 0)
	ov2:=mol.Coord(orient_atoms[1], 0)
	//now we center the thing in the beta carbon of D124
	err=SubRow(mol.Coords[0],ov2)
	//Now the rotation
	ov1=mol.Coord(orient_atoms[0], 0) //make sure we have the correct versions
	ov2=mol.Coord(orient_atoms[1], 0)  //same
	orient:=ov2.Copy()	
	orient.SubtractDense(ov1)
	rotation,err:=GetSwitchZ(mol.Coords[0],orient)
	fmt.Println("rotation: ",rotation)
	if err!= nil {
		Te.Error(err)
		}
	mol.Coords[0]=matrix.ParallelProduct(mol.Coords[0],rotation)
	fmt.Println(orient_atoms[1], mol.Atoms[orient_atoms[1]],mol.Atoms[orient_atoms[0]])//, mol.Coords[0][orient_atoms[1]])
	if err!=nil{
		Te.Error(err)
		}
	PdbWrite(&mol,"test/2c9v-aligned.pdb")
	}


//TestGeo opens the sample.xyz file in the test directory, and pull a number of hardcoded atoms
//In the direction of a hardcoded vectos. It builds 12 files with the pulled atoms  displaced by
//different along the pulling vector
func TestGeo(Te *testing.T) {
	pulled_atoms:=[7]int{43,41,42,40,85,86,87}
	pulling_vector:=[2]int{40,88}
	var mol Molecule
	a,b,err:=XyzRead("test/sample.xyz")
	if err!=nil{
		Te.Error(err)
		}
	mol.Atoms=a
	mol.Coords=b
	pulled_res,err:=mol.GetCoords(pulled_atoms[:], 0)
	if err!=nil{
		Te.Error(err)
		}
	at1:=mol.Coord(pulling_vector[0],0)
	vector:=mol.Coord(pulling_vector[1],0)
	vector=vector.Copy()
	err=vector.SubtractDense(at1)
	if err!=nil{
		Te.Error(err)
		}
	vector=Unitarize(vector)
	var scale_factors = [12]float64{-1.0, -2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}
	for _,scaling:=range(scale_factors){
		vec:=vector.Copy()
		pulled:=pulled_res.Copy()
		vec.Scale(scaling)
		err=AddRow(pulled,vec)
		if err!=nil{
			Te.Error(err)
			}
		mol.SetCoords(pulled_atoms[:], 0, pulled)
		err=mol.Corrupted()
		if err!=nil{
			Te.Error(err)
			}
		XyzWrite(&mol, 0, strings.Replace("test/sample_xxxx.xyz","xxxx",strconv.FormatFloat(scaling, 'f', 1, 64),1)) //There might be an easier way of creating the filenames
	}
	//fmt.Println(mol.Atoms,mol.Coords,err,pulled,vector,vector.TwoNorm())
	}
