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
//import "github.com/skelterjohn/go.matrix"
//import "fmt"
import "testing"

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
	pulled,err:=mol.GetCoords(pulled_atoms[:], 0)
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
	vector.Scale(3)
	err=AddRow(pulled,vector)
	if err!=nil{
		Te.Error(err)
		}
	mol.SetCoords(pulled_atoms[:], 0, pulled)
	err=mol.Corrupted()
	if err!=nil{
		Te.Error(err)
		}
	XyzWrite(&mol, 0, "test/sample_mod.xyz") 
	
	//fmt.Println(mol.Atoms,mol.Coords,err,pulled,vector,vector.TwoNorm())
	}
