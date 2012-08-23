/*
 * general_test.go
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */



package chem
import "github.com/skelterjohn/go.matrix"
import "fmt"
import "testing"
import "strings"
import "strconv"

func Distance(first, second *matrix.DenseMatrix)(float64, error){
	if first.Rows() != 1 || second.Rows() != 1 || first.Cols() != 3 || second.Cols()!=3{
		return 0.0, fmt.Errorf("Ill-formated input matrices")
		} 
	distance:=second.Copy()
	_=SubRow(distance,first) //We already checked the format so we shouldnt need to check errors here
	dist:=ditance.TwoNorm()
	return dist, nil
	}

func TestFilterWaters(Te *testing.T) {
	var mol Molecule
	ats,coords,bfac,err:=PdbRead("test/2c9v.pdb",true)
	if err!=nil{
		Te.Error(err)
		}
	mol.Atoms=ats
	mol.Coords=coords
	mol.Bfactors=bfac
 	//Now comes the fun
	frames:=len(mol.Coords)
	err:=mol.Corrupted
	if err!=nil{
		Te.Error(err)
		}
	atoms_per_frame:=len(mol.Coords[0])
	var coord_saved:= make([]int)
	var oxygens:= make([]int,0,2)
	d124:=false //if we had the molecule already
	site:=false //this is the last molecule of the active site
	for j,k := range(mol.Atoms){
		if mol.Atoms[j].Molid==124 && d124==false{
			coords_saved=append(coords_saved,j)
			if mol.Atoms[j].Name="OD1" || mol.Atoms[j].Name="OD2"{
				oxygen=append(oxygen,j)
				}
			}
		if mol.Atoms[j].Molid==125 && d124==false{
			d124=true
			}
		if (mol.Atoms[j].Molid==46 || mol.Atoms[j].Molid==48 || mol.Atoms[j].Molid==63 || mol.Atoms[j].Molid==71 || mol.Atoms[j].Molid==80 || mol.Atoms[j].Molid==83 || mol.Atoms[j].Molid==120 || mol.Atoms[j].Name=="ZN" || mol.Atoms[j].Name=="CU") && site==false{
			coords_saved=append(coords_saved,j)
			}
		if mol.Atoms[j].Molid==156{ // I might need to change this number	
			site=true //This mean that we have already collected the active site atoms
			}	
		}
	//Now we look for the 8 waters closer to the ASP124 oxygens
	var waters molecule
	for i:=0;i<frames;i++{
		waters.Coords=append(waters.Coords,matrix.MakeDenseMatrix) //Should work
		waters.Bfactors=append(waters.Bfactors,make([]float64)) //Should work
		latest:=len(waters.Coords)-1
		od1:=mol.Coords.GetRowVecto(oxygen[0])
		od2:=mol.Coords.GetRowVecto(oxygen[1])
		toinclude:=8 //number of water molecules to be included
		for j:=0,j<mol.Coords[i].Rows();j++{
			if mol.Atoms[j].Molname=="SOL"{
				if i==0{
					waters.Atoms= append(water.Atoms,mol.Atoms[j])
					}
				waters.Coords[latest]=append(water.Coords[latest],mol.Coords[i][j])
				waters.Bfactors[latest]=append(water.Bfactors[latest],mol.Bfactors[i][j])
				}
			}
		}
	
	}

//TestChangeAxis
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
