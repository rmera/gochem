/*
 * json.go, part of gochem.
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

import (
	"bufio"
	"encoding/json"
	"strings"
	"fmt"
	"os"
)

type JSONAtom struct {
	A      *Atom
	Coords []float64
	Bfac   float64
}

type JSONCoords struct {
	Coords []float64
}

//An easily JSON-serializable error type, 
type JSONError struct {
	IsError      bool //If this is false (no error) all the other fields will be at their zero-values.
	InOptions    bool //If error, was it in parsing the options?
	InSelections bool //Was it in parsing selections?
	InProcess    bool
	InPostProcess bool //was it in preparing the output?
	Selection    string //Which selection?
	State        int    //Which state of it?
	Atom         int
	Function     string //which go function gave the error
	Message      string //the error itself
}
//implements the error interface
func (J *JSONError) Error () string{
	return  J.Message
}

//Serializes the error. Panics on failure.
func (J *JSONError) Marshal() []byte{
	ret, err2 := json.Marshal(J)
	if err2 != nil {
		panic(strings.Join([]string{J.Error(), err2.Error()}, " - ")) //well, shit.
	}
	return ret
}

//Information to be passed back to the calling program.
type JSONInfo struct {
	Molecules int
	Bfactors  bool
	SS        bool
	FramesPerMolecule []int
	AtomsPerMolecule  []int
	FloatInfo  [][]float64
	StringInfo [][]string
	IntInfo    [][]int
}
func (J *JSONInfo) Marshal() ([]byte,*JSONError){
	ret,err:=json.Marshal(J)
	if err!=nil{
		return nil, MakeJSONError("postprocess","JSONInfo.Marshal",err)
	}
	return ret, nil
}



//Options passed from the calling external program
type JSONOptions struct {
	SelNames      []string
	AtomsPerSel   []int //Knowing in advance makes memory allocation more efficient
	StatesPerSel  []int //How many snapshots a traj has?
	StringOptions [][]string
	IntOptions    [][]int
	BoolOptions   [][]bool
	FloatOptions  [][]float64
}



//Takes an error and some additional info to create a JSON error
func MakeJSONError(where, function string, err error) *JSONError {
	jerr := new(JSONError)
	jerr.IsError = true
	switch where {
	case "options":
		jerr.InOptions = true
	case "selection":
		jerr.InSelections = true
	case "postprocess":
		jerr.InPostProcess = true
	default:
		jerr.InProcess = true
	}
	jerr.Function = function
	jerr.Message = err.Error()
	return jerr
}



func DecodeJSONOptions(stdin *bufio.Reader) (*JSONOptions, *JSONError) {
	line, err := stdin.ReadBytes('\n')
	if err != nil {
		return nil, MakeJSONError("options","DecodeJSONOptions",err)
	}
	ret := new(JSONOptions)
	err = json.Unmarshal(line, ret)
	if err != nil {
		return nil, MakeJSONError("options","DecodeJSONOptions",err)
	}
	return ret, nil
}


//Decodes a JSON molecule into a gochem molecule. Can handle several frames (all of which need to have the same amount of atoms). It does
//not collect the b-factors.
func DecodeJSONMolecule(stream *bufio.Reader, atomnumber, frames int) (*Topology, []*CoordMatrix, *JSONError) {
	atoms := make([]*Atom, 0, atomnumber)
	coordset := make([]*CoordMatrix, 0, frames)
	rawcoords := make([]float64, 0, 3*atomnumber)
	for i := 0; i < atomnumber; i++ {
		line, err := stream.ReadBytes('\n') //Using this function allocates a lot without need. There is no function that takes a []bytes AND a limit. I might write one at some point.
		if err != nil {
			break
		}
		at := new(Atom)
		err = json.Unmarshal(line, at)
		if err != nil {
			return nil, nil, MakeJSONError("selection","DecodeJSONMolecule",err)
		}
		atoms = append(atoms, at)
		line, err = stream.ReadBytes('\n') //See previous comment.
		if err != nil {
			break
		}
		ctemp := new(JSONCoords)
		if err = json.Unmarshal(line, ctemp); err != nil {
			return nil, nil, MakeJSONError("selection","DecodeJSONMolecule",err)
		}
		rawcoords = append(rawcoords, ctemp.Coords...)
	}
	mol, _ := MakeTopology(atoms, 0, 0) //no idea of the charge or multiplicity
	coords := NewCoords(rawcoords)
	coordset = append(coordset, coords)
	if frames == 1 {
		return mol, coordset, nil
	}
	for i:=0;i<(frames-1);i++{
		coords,err:=DecodeJSONCoords(stream,atomnumber)
		if err!=nil{
			return mol, coordset, MakeJSONError("selection","DecodeJSONMolecule",fmt.Errorf("Error reading the %d th frame: %s",i+2,err.Error()))
			}
		coordset=append(coordset,coords)
	}
	return mol, coordset, nil

}

func DecodeJSONCoords(stream *bufio.Reader, atomnumber int) (*CoordMatrix, *JSONError) {
	rawcoords := make([]float64, 0, 3*atomnumber)
	for i := 0; i < atomnumber; i++ {
		line, err := stream.ReadBytes('\n') //Using this function allocates a lot without need. There is no function that takes a []bytes AND a limit. I might write one at some point.
		if err != nil {
			break
		}
		ctemp := new(JSONCoords)
		if err = json.Unmarshal(line, ctemp); err != nil {
			return nil, MakeJSONError("selection","DecodeJSONCoords",err)
		}
		rawcoords = append(rawcoords, ctemp.Coords...)
	}
	coords := NewCoords(rawcoords)
	return coords, nil
}




func TransmitMoleculeJSON(mol Atomer, coordset, bfactors []*CoordMatrix, ss [][]string, out *os.File) *JSONError{
	jmol,err:=EncodeAtoms2JSON(mol)
	if err!=nil{
		return err
	}
	for _,i:=range(jmol){
		fmt.Fprint(out,i)
		fmt.Fprint(out,"\n")
	}
	for _,coords:=range(coordset){
		jcoords,err:=EncodeCoords2JSON(coords)
		if err!=nil{
			return err
		}
		for _,i:=range(jcoords){
			fmt.Fprint(out,i)
			fmt.Fprint(out,"\n")
		}
	}
	if bfactors!=nil{
		jb:=new(jSONbfac)
		for _,b:=range(bfactors){
			jb.Bfactors=b.Col(nil,0)
			jjb,err2:=json.Marshal(jb)
			if err2!=nil{
				return MakeJSONError("postprocess","TransmitMoleculeJson(bfactors)",err2)
			}
			fmt.Fprint(out,jjb)
			fmt.Fprint(out,"\n")
		}
	}
	if ss!=nil{
		jss:=new(jSONss)
		for _,s:=range(ss){
			jss.SS=s
			jjss,err2:=json.Marshal(jss)
			if err2!=nil{
				return MakeJSONError("postprocess","TransmitMoleculeJson(ss)",err2)
			}
			fmt.Fprint(out,jjss)
			fmt.Fprint(out,"\n")
		}
	}
	return nil
}

type jSONbfac struct {
	Bfactors []float64
}

type jSONss struct {
	SS []string
}


type jSONCoords struct{
	coords []float64
}

func EncodeAtoms2JSON(mol Atomer) ([][]byte, *JSONError){
	ret:=make([][]byte,0,mol.Len())
	for i:=0;i<mol.Len();i++{
		b,err:=json.Marshal(mol.Atom(i))
		if err!=nil{
			return  nil, MakeJSONError("postprocess","EncodeAtoms2JSON",err)
		}
		ret=append(ret,b)
	}
	return ret, nil
}



func EncodeCoords2JSON(coords *CoordMatrix)  ([][]byte, *JSONError){
	c:=new(jSONCoords)
	t:=make([]float64,3,3)
	ret:=make([][]byte,0,coords.NumVec())
	for i:=0;i<coords.NumVec();i++{
		c.coords=coords.Row(t,i)
		b,err:=json.Marshal(c)
		if err!=nil{
			return  nil, MakeJSONError("postprocess","EncodeCoords2JSON",err)
		}
		ret=append(ret,b)
	}
	return ret, nil
}



