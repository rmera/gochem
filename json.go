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
	"os"
)

type pymolatom struct {
	A      *Atom
	Coords []float64
	Bfac   float64
}

type JsonCoords struct {
	Coords []float64
}

func EncodeAtoms(mol Atomer, coord, bfact *CoordMatrix) error {
	enc := json.NewEncoder(os.Stdout)

	for i := 0; i < mol.Len(); i++ {
		//		err := enc.Encode(coord.Row(i))
		//		if err != nil {
		//			return err
		//		}
		err := enc.Encode(pymolatom{mol.Atom(i), coord.Row(i), bfact.At(i, 0)})
		if err != nil {
			return err
		}
	}
	return nil
}

func EncodeCoords(mol Atomer, coord, bfact *CoordMatrix) error {
	enc := json.NewEncoder(os.Stdout)

	for i := 0; i < mol.Len(); i++ {
		//		err := enc.Encode(coord.Row(i))
		//		if err != nil {
		//			return err
		//		}
		err := enc.Encode(coord.Row(i))
		if err != nil {
			return err
		}
	}
	return nil
}

func DecodePyMolStream(stream *bufio.Reader) (*Topology, *CoordMatrix, error) {
	var err error
	atoms := make([]*Atom, 0, 10)
	rawcoords := make([]float64, 0, 30)
	for i := 0; err == nil; i++ {
		line, err := stream.ReadBytes('\n')
		if err != nil {
			break
		}
		at := new(Atom)
		err = json.Unmarshal(line, at)
		if err != nil {
			return nil, nil, err
		}
		atoms = append(atoms, at)
		line, err = stream.ReadBytes('\n')
		if err != nil {
			break
		}
		ctemp := new(JsonCoords)
		if err = json.Unmarshal(line, ctemp); err != nil {
			return nil, nil, err
		}
		rawcoords = append(rawcoords, ctemp.Coords...)
	}
	mol, _ := MakeTopology(atoms, 0, 0) //no idea of the charge or multiplicity
	coords := NewCoords(rawcoords)
	return mol, coords, nil
}
