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
//	"os"
)

type JSONAtom struct {
	A      *Atom
	Coords []float64
	Bfac   float64
}

type JSONCoords struct {
	Coords []float64
}

type JSONError struct {
	Error    bool //If this is false (no error) all the other fields will be at their zero-values.
	InOptions bool //If error, was it in parsing the options?
	InSelections bool //Was it in parsing selections?
	InProcess bool
	Selection  string  //Which selection?
	State   int    //Which state of it?
	Atom  int     
	Function string  //which go function gave the error
	Message string   //the error itself
}


type JSONOptions struct {
	SelNames []string
	AtomsPerSel []int  //Knowing in advance makes memory allocation more efficient
	StatesPerSel []int   //How many snapshots a traj has?
	StringOptions [][]string
	IntOptions    [][]int
	BoolOptions   [][]bool
	FloatOptions  [][]float64
}


func MakeJSONError(where, function string, err error)([]byte){
	jerr:=new(JSONError)
	jerr.Error=true
	switch where{
	case "options":
		jerr.InOptions=true
	case "selection":
		jerr.InSelections=true
	default:
		jerr.InProcess=true
	}
	jerr.Function=function
	jerr.Message=err.Error()
	ret,err2:=json.Marshal(jerr)
	if err2!=nil{
		panic(strings.Join([]string{err.Error(),err2.Error()}," - ")) //well, shit.
		}
	return ret
}


func DecodeJSONOptions(stdin *bufio.Reader) (*JSONOptions, error){
	line,err:=stdin.ReadBytes('\n')
	if err!=nil{
		return nil, err
	}
	ret:=new(JSONOptions)
	err=json.Unmarshal(line, ret)
	if err!=nil{
		return nil, err
	}
	return ret, nil
}


func DecodeJSONMolecule(stream *bufio.Reader, atomnumber int) (*Topology, *CoordMatrix, error) {
	atoms := make([]*Atom, 0, atomnumber)
	rawcoords := make([]float64, 0, 3*atomnumber)
	for i := 0; i<atomnumber; i++ {
		line, err := stream.ReadBytes('\n')  //Using this function allocates a lot without need. There is no function that takes a []bytes AND a limit. I might write one at some point.
		if err != nil {
			break
		}
		at := new(Atom)
		err = json.Unmarshal(line, at)
		if err != nil {
			return nil, nil, err
		}
		atoms = append(atoms, at)
		line, err = stream.ReadBytes('\n')   //See previous comment.
		if err != nil {
			break
		}
		ctemp := new(JSONCoords)
		if err = json.Unmarshal(line, ctemp); err != nil {
			return nil, nil, err
		}
		rawcoords = append(rawcoords, ctemp.Coords...)
	}
	mol, _ := MakeTopology(atoms, 0, 0) //no idea of the charge or multiplicity
	coords := NewCoords(rawcoords)
	return mol, coords, nil
}


func DecodeJSONCoords(stream *bufio.Reader, atomnumber int) (*CoordMatrix, error) {
	rawcoords := make([]float64, 0, 3*atomnumber)
	for i := 0; i<atomnumber; i++ {
		line, err := stream.ReadBytes('\n')  //Using this function allocates a lot without need. There is no function that takes a []bytes AND a limit. I might write one at some point.
		if err != nil {
			break
		}
		ctemp := new(JSONCoords)
		if err = json.Unmarshal(line, ctemp); err != nil {
			return nil, err
		}
		rawcoords = append(rawcoords, ctemp.Coords...)
	}
	coords := NewCoords(rawcoords)
	return coords, nil
}

