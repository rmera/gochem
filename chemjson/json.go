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

package chemjson

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

//A ready-to-serialize container for an atom.
type Atom struct {
	A      *chem.Atom
	Coords []float64
	Bfac   float64
}

//A ready-to-serialize container for coordinates
type Coords struct {
	Coords []float64
}

//An easily JSON-serializable error type,
type Error struct {
	deco          []string
	IsError       bool //If this is false (no error) all the other fields will be at their zero-values.
	InOptions     bool //If error, was it in parsing the options?
	InSelections  bool //Was it in parsing selections?
	InProcess     bool
	InPostProcess bool   //was it in preparing the output?
	Selection     string //Which selection?
	State         int    //Which state of it?
	Atom          int
	Function      string //which go function gave the error
	Message       string //the error itself
}

//Error implements the error interface
func (J *Error) Error() string {
	return J.Message
}

//Decorate will add the dec string to the decoration slice of strings of the error,
//and return the resulting slice.
func (err Error) Decorate(dec string) []string {
	if dec == "" {
		return err.deco
	}
	err.deco = append(err.deco, dec)
	return err.deco
}

//Serializes the error. Panics on failure.
func (J *Error) Marshal() []byte {
	ret, err2 := json.Marshal(J)
	if err2 != nil {
		panic(strings.Join([]string{J.Error(), err2.Error()}, " - ")) // Yo, dawg, I heard you like errors, so I got an error while serializing your error so you can... you know the drill.
	}
	return ret
}

//Information to be passed back to the calling program.
type Info struct {
	Molecules         int
	Bfactors          bool
	SS                bool
	FramesPerMolecule []int
	AtomsPerMolecule  []int
	FloatInfo         [][]float64
	StringInfo        [][]string
	IntInfo           [][]int
	BoolInfo          [][]bool
	Energies          []float64
}

//Send Marshals the info and writes to out, returns an error or nil
func (J *Info) Send(out io.Writer) *Error {
	enc := json.NewEncoder(out)
	if err := enc.Encode(J); err != nil {
		return NewError("postprocess", "JSONInfo.Marshal", err)
	}
	return nil
}

//Options passed from the calling external program
type Options struct {
	SelNames      []string
	AtomsPerSel   []int //Knowing in advance makes memory allocation more efficient
	StatesPerSel  []int //How many snapshots a traj has?
	StringOptions [][]string
	IntOptions    [][]int
	BoolOptions   [][]bool
	FloatOptions  [][]float64
}

//Takes an error and some additional info to create a json-marshal-ble error
func NewError(where, function string, err error) *Error {
	jerr := new(Error)
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

//DecodeOptions Decodes or unmarshals json options into an Options structure
func DecodeOptions(stdin *bufio.Reader) (*Options, *Error) {
	line, err := stdin.ReadBytes('\n')
	if err != nil {
		return nil, NewError("options", "DecodeOptions", err)
	}
	ret := new(Options)
	err = json.Unmarshal(line, ret)
	if err != nil {
		return nil, NewError("options", "DecodeOptions", err)
	}
	return ret, nil
}

//DecodeMolecule Decodes a JSON molecule into a gochem molecule. Can handle several frames (all of which need to have the same amount of atoms). It does
//not collect the b-factors.
func DecodeMolecule(stream *bufio.Reader, atomnumber, frames int) (*chem.Topology, []*v3.Matrix, *Error) {
	const funcname = "DecodeMolecule" //for the error
	atoms := make([]*chem.Atom, 0, atomnumber)
	coordset := make([]*v3.Matrix, 0, frames)
	rawcoords := make([]float64, 0, 3*atomnumber)
	for i := 0; i < atomnumber; i++ {
		line, err := stream.ReadBytes('\n') //Using this function allocates a lot without need. There is no function that takes a []bytes AND a limit. I might write one at some point.
		if err != nil {
			break
		}
		at := new(chem.Atom)
		err = json.Unmarshal(line, at)
		if err != nil {
			return nil, nil, NewError("selection", funcname, err)
		}
		atoms = append(atoms, at)
		line, err = stream.ReadBytes('\n') //See previous comment.
		if err != nil {
			break
		}
		ctemp := new(Coords)
		if err = json.Unmarshal(line, ctemp); err != nil {
			return nil, nil, NewError("selection", funcname, err)
		}
		rawcoords = append(rawcoords, ctemp.Coords...)
	}
	mol := chem.NewTopology(-1, 99999, atoms) //no idea of the charge or multiplicity
	coords, err := v3.NewMatrix(rawcoords)
	if err != nil {
		return nil, nil, NewError("selection", funcname, err)
	}
	coordset = append(coordset, coords)
	if frames == 1 {
		return mol, coordset, nil
	}
	for i := 0; i < (frames - 1); i++ {
		coords, err := DecodeCoords(stream, atomnumber)
		if err != nil {
			return mol, coordset, NewError("selection", funcname, fmt.Errorf("Error reading the %d th frame: %s", i+2, err.Error()))
		}
		coordset = append(coordset, coords)
	}
	return mol, coordset, nil

}

//Decodecoords decodes streams from a bufio.Reader containing 3*atomnumber JSON floats into a v3.Matrix with atomnumber rows.
func DecodeCoords(stream *bufio.Reader, atomnumber int) (*v3.Matrix, *Error) {
	const funcname = "DecodeCoords"
	rawcoords := make([]float64, 0, 3*atomnumber)
	for i := 0; i < atomnumber; i++ {
		line, err := stream.ReadBytes('\n') //Using this function allocates a lot without need. There is no function that takes a []bytes AND a limit. I might write one at some point.
		if err != nil {
			break
		}
		ctemp := new(Coords)
		if err = json.Unmarshal(line, ctemp); err != nil {
			return nil, NewError("selection", funcname, err)
		}
		rawcoords = append(rawcoords, ctemp.Coords...)
	}
	coords, err := v3.NewMatrix(rawcoords)
	if err != nil {
		return nil, NewError("selection", funcname, err)
	}
	return coords, nil
}

//Takes a chem.Atomer, coordinates, bfactos and secondary strucuctures, encodes them and writes them to the given io.writer
func SendMolecule(mol chem.Atomer, coordset, bfactors []*v3.Matrix, ss [][]string, out io.Writer) *Error {
	const funcname = "SendMolecule"
	enc := json.NewEncoder(out)
	if err := EncodeAtoms(mol, enc); err != nil {
		return err
	}
	for _, coords := range coordset {
		if err := EncodeCoords(coords, enc); err != nil {
			return err
		}
	}
	if bfactors != nil {
		jb := new(jSONbfac)
		for _, b := range bfactors {
			jb.Bfactors = b.Col(nil, 0)
			if err := enc.Encode(jb); err != nil {
				return NewError("postprocess", funcname+"(bfactors)", err)
			}
		}
	}
	if ss != nil {
		jss := new(jSONss)
		for _, s := range ss {
			jss.SS = s
			if err := enc.Encode(jss); err != nil {
				return NewError("postprocess", funcname+"(ss)", err)
			}
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

type jSONCoords struct {
	Coords []float64
}

//Encodes a goChem Atomer into a JSON
func EncodeAtoms(mol chem.Atomer, enc *json.Encoder) *Error {
	const funcname = "EncodeAtoms"
	if mol == nil {
		return nil //Its assumed to be intentional.
	}
	for i := 0; i < mol.Len(); i++ {
		if err := enc.Encode(mol.Atom(i)); err != nil {
			return NewError("postprocess", funcname, err)
		}
	}
	return nil
}

//Encodes a set of coordinates into JSON
func EncodeCoords(coords *v3.Matrix, enc *json.Encoder) *Error {
	c := new(Coords)
	t := make([]float64, 3, 3)
	for i := 0; i < coords.NVecs(); i++ {
		c.Coords = coords.Row(t, i)
		if err := enc.Encode(c); err != nil {
			return NewError("postprocess", "chemjson.EncodeCoords", err)
		}
	}
	return nil
}
