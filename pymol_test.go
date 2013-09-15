// +build pymol 

/*
 * pymol_test.go
 *
 * Copyright 2013  <rmera@Holmes>
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
import (
	"encoding/json"
	"fmt"
	"testing"
)

func TestAtomEncode(Te *testing.T) {
	mol, err := PDBRead("test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	err = EncodeAtoms(mol, mol.Coords[0], mol.Bfactors[0])
	if err != nil {
		Te.Error(err)
	}
	a := []byte(`{"Name": "ZN", "Chain": "F", "SS": "", "Symbol": "ZN", "Molname": "ZN", "Id": 2427, "Molid": 155}`)
	at := new(Atom)
	json.Unmarshal(a, at)
	fmt.Println(at)
	//for the 3 residue  I should get -131.99, 152.49.
}
