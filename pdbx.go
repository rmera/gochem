/*
 * mmcif.go, part of gochem.
 *
 *
 * Copyright 2024 rmeraaatacademicosdotutadotcl
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
 * goChem is developed at Universidad de Tarapaca (UTA)
 *
 *
 */

package chem

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	v3 "github.com/rmera/gochem/v3"
)

var tl func(string) string = strings.ToLower

// PDBRRead reads a pdb file from an io.Reader. Returns a Molecule. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not
// really well set up right now.
// read_additional is now "deprecated", it will be set to true regardless. I have made it into
func PDBxRead(pdb io.Reader) (*Molecule, error) {
	bufiopdb := bufio.NewReader(pdb)
	mol, err := pdbBufIORead(bufiopdb)
	return mol, errDecorate(err, "PDBReaderREad")
}

// PDBFileRead reads a pdb file from an io.Reader. Returns a Molecule. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not
// really well set up right now. read_additional is now deprecated. The reader will just read
func PDBxFileRead(pdbname string) (*Molecule, error) {
	pdbxfile, err := os.Open(pdbname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, err
	}
	defer pdbxfile.Close()
	pdb := bufio.NewReader(pdbxfile)
	mol, err := pdbxBufIORead(pdb)
	return mol, err
}

func pdbxNextLoop(pdb *bufio.Reader) (*bufio.Reader, string, error) {
	for {
		line, err := pdb.ReadString('\n')
		if err != nil {
			return pdb, line, err
		}
		if strings.HasPrefix(tl(line), "loop_") {
			return pdb, line, nil
		}
	}

}

type pdbxmap map[string]int

// adds i to the map[string] entry, if it exists. If not,
// does nothing. Returns the map.
func (m pdbxmap) add(s string, i int) pdbxmap {
	s = strings.TrimSpace(s)
	s = strings.Replace(s, "\n", "", -1)
	if _, ok := m[s]; ok {
		m[s] = i
	}
	return m
}

// returns the integaer corresponding to the given string in the map
// or -1 if the string is not a key in the map.
func (m pdbxmap) get(s string) int {
	if i, ok := m[s]; ok {
		return i
	}
	return -1
}

func pdbxFillAtom(at *Atom, data []string, m pdbxmap) error {
	do := func(at *Atom, s string, data []string, m pdbxmap, f func(*Atom, string) error) error {
		s = strings.TrimSpace(s) //just in case I make a mistake
		//	fmt.Println(data, m, s)
		k := m.get(s)
		if k >= 0 {
			if k >= len(data) {
				return fmt.Errorf("Index out of range: %d, %v", k, data)
			}
			//	println("DATA", k, data[k]) /////////////////////
			err := f(at, data[k])
			return err
		}
		return nil
	}
	//We start with the simpler string fields
	//symbol
	f := func(a *Atom, s string) error {
		a.Symbol = s
		return nil
	}
	do(at, "_atom_site.type_symbol", data, m, f)
	//name
	f = func(a *Atom, s string) error {
		a.Name = s
		return nil
	}
	//println("NAME", at.Name) /////////////////////////
	do(at, "_atom_site.auth_atom_id", data, m, f)
	if at.Symbol == "" {
		at.Symbol, _ = symbolFromName(at.Name)
	}
	//molname
	f = func(a *Atom, s string) error {
		a.MolName = s
		return nil
	}
	do(at, "_atom_site.auth_comp_id", data, m, f)
	at.MolName1 = three2OneLetter[at.MolName]

	//char16
	f = func(a *Atom, s string) error {
		a.Char16 = s[0]
		return nil
	}
	do(at, "_atom_site.label_alt_id", data, m, f)
	//chain
	f = func(a *Atom, s string) error {
		a.Chain = s
		return nil
	}
	do(at, "_atom_site.auth_asym_id", data, m, f)

	//Now the integer fields
	//ID
	f = func(a *Atom, s string) error {
		n, err := strconv.Atoi(s)
		if err != nil {
			return fmt.Errorf("pdbxFillAtom: Couldn't parse ID from %s: %w", s, err)
		}
		a.ID = n
		return nil
	}
	err := do(at, "_atom_site.id", data, m, f)
	if err != nil {
		return err
	}
	//molID
	f = func(a *Atom, s string) error {
		n, err := strconv.Atoi(s)
		if err != nil {
			return fmt.Errorf("pdbxFillAtom: Couldn't parse MolID from %s: %w", s, err)
		}
		a.MolID = n
		return nil
	}
	err = do(at, "_atom_site.auth_seq_id", data, m, f)
	if err != nil {
		return err
	}
	//Now the floating point fields, except for the coordinates and b-factors
	//occupancy
	f = func(a *Atom, s string) error {
		n, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return fmt.Errorf("pdbxFillAtom: Couldn't parse Occupancy from %s: %w", s, err)
		}
		a.Occupancy = n
		return nil
	}
	err = do(at, "_atom_site.occupancy", data, m, f)
	if err != nil {
		return err
	}
	//Charge, but we won't do anything if we somehow can't read it.
	f = func(a *Atom, s string) error {
		n, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return nil
		}
		a.Charge = n
		return nil
	}
	do(at, "_atom_site.pdbx_formal_charge ", data, m, f)
	//And, finally, the boolean field
	f = func(a *Atom, s string) error {
		if s != "ATOM" {
			a.Het = true
		} else {
			a.Het = false
		}
		return nil
	}
	do(at, "_atom_site.group_pdb", data, m, f)
	return nil
}

func pdbxFillBfac(data []string, bf []float64, m pdbxmap) ([]float64, error) {
	v := "_atom_site.b_iso_or_equiv"
	if m[v] >= 0 && m[v] < len(data) {
		fl, err := strconv.ParseFloat(data[m[v]], 64)
		if err != nil {
			return bf, fmt.Errorf("pdbxFillBfac: Couldn't parse bfactor from %s: %w", data[m[v]], err)
		}
		bf = append(bf, fl)
	} else {
		return bf, fmt.Errorf("pdbxFillBfac: Field %s not present in data %v", v, data)
	}
	return bf, nil
}

func pdbxFillCoords(data []string, coord []float64, m pdbxmap) ([]float64, error) {
	c := []string{"_atom_site.cartn_x", "_atom_site.cartn_y", "_atom_site.cartn_z"}
	for j, v := range c {
		if m[v] >= 0 && m[v] < len(data) {
			fl, err := strconv.ParseFloat(data[m[v]], 64)
			if err != nil {
				return coord, fmt.Errorf("pdbxFillCoord: Couldn't parse %d cartesian coordinate from %s: %w", j, data[m[v]], err)
			}
			coord = append(coord, fl)
		} else {
			return coord, fmt.Errorf("pdbxFillCoord: Field %s not present in data %v", v, data)
		}
	}
	return coord, nil
}

func pdbxBufIORead(pdb *bufio.Reader) (*Molecule, error) {
	m := pdbxmap(ma)
	molecule := make([]*Atom, 0)
	coords := make([][]float64, 1, 1)
	coords[0] = make([]float64, 0, 3)
	bfactors := make([][]float64, 1, 1)
	bfactors[0] = make([]float64, 0)
	currentmodel := 1
	var reading bool
	var field int = 0
	var havebfactors bool = true
	hp := strings.HasPrefix
	trimall := func(s string) string { return strings.TrimSpace(strings.Replace(s, "\n", "", -1)) }

	for {
		line, err := pdb.ReadString('\n')
		if err != nil {
			break
		}
		if hp(line, "#") || hp(line, ";") || trimall(line) == "" {
			continue
		}
		if !reading && strings.HasPrefix(tl(line), "_atom_site") {
			reading = true
			field = 0
		}
		if !reading {
			pdb, line, err = pdbxNextLoop(pdb)
			if err != nil {
				break
			}
			continue
		}
		if strings.HasPrefix(line, "loop_") { //new section
			reading = false
			continue
		}
		// We shouldn't be here if reading is false
		if strings.HasPrefix(line, "_") {
			if !hp(line, "_atom_site") || hp(line, "_atom_site_anisotrop") { //a new section started
				reading = false
				continue
			}
			//	fmt.Println("This is reading!!", line, field) ///////////

			m.add(tl(line), field)
			field++
		} else {
			//Here we should be reading the content lines.
			//we first see if we have a model number, and whether
			fields := strings.Fields(line)
			modkey := m.get("_atom_site.pdbx_pdb_model_num")
			if modkey >= 0 {
				if modkey >= len(fields) {
					return nil, fmt.Errorf("pdbxBufIORead: Model field out of range: %d (index) %d (len),  %v", modkey, len(fields), fields)
				}
				mod := fields[modkey]
				model, err := strconv.Atoi(mod)
				if err != nil {
					return nil, fmt.Errorf("pdbxBufIORead: Couldn't parse model number from %s: %w", mod, err)
				}
				if model > currentmodel {
					nats := len(coords[len(coords)-1])
					coords = append(coords, make([]float64, 0, nats))
					bfactors = append(bfactors, make([]float64, 0, nats))
					currentmodel = model
				}
			}
			//we don't read the atoms again for the next models.
			if currentmodel == 1 {
				at := new(Atom)
				err := pdbxFillAtom(at, fields, m)
				if err != nil {
					return nil, fmt.Errorf("pdbxBufIORead: Couldn't read atom %d: %w", len(molecule)+1, err)
				}
				molecule = append(molecule, at)
			}
			//The following we always try to read.
			c := len(coords) - 1
			b := len(bfactors) - 1 //c and b should be equal, but not sure if its critical enough to check it.
			coords[c], err = pdbxFillCoords(fields, coords[c], m)
			if err != nil {
				return nil, fmt.Errorf("pdbxBufIORead: Couldn't read %d th coordinates for frame %d: %w", len(coords[c])+1, currentmodel, err)
			}
			bfactors[b], err = pdbxFillBfac(fields, bfactors[b], m)
			if err != nil && havebfactors {
				//I might remove this, but not before more testing.
				//It can very well be that the PDB just doesn't contain b-factors.
				log.Printf("pdbxBufIORead: Couldn't read %d th bfactors for frame %d: %v", len(bfactors[b])+1, currentmodel, err)
				havebfactors = false

			}
		}

	}
	top := NewTopology(0, 1, molecule)
	var err error
	frames := len(coords)
	mcoords := make([]*v3.Matrix, frames, frames) //Final thing to return
	for i := 0; i < frames; i++ {
		mcoords[i], err = v3.NewMatrix(coords[i])
		if err != nil {
			return nil, fmt.Errorf("pdbBufIORead: Couldn't transfor coordinates from frame %d: %w", i, err)
		}
	}
	//if something happened during the process
	if err != nil {
		return nil, fmt.Errorf("pdbBufIORead: %w", err)

	}
	if !havebfactors {
		bfactors = nil
	}
	returned, err := NewMolecule(mcoords, top, bfactors)
	if err != nil {
		return returned, fmt.Errorf("pdbBufIORead: %w", err)
	}
	return returned, nil
}

func PDBxFileWrite(name string, coords []*v3.Matrix, mol Atomer, bfact [][]float64) error {
	pdb, err := os.Create(name)
	if err != nil {
		return fmt.Errorf("PDBxFileWrite: %w", err)
	}
	defer pdb.Close()
	return PDBxWrite(pdb, coords, mol, bfact, false, strings.Replace(name, ".pdb", "", -1))
}

func PDBxCompactFileWrite(name string, coords []*v3.Matrix, mol Atomer, bfact [][]float64) error {
	pdb, err := os.Create(name)
	if err != nil {
		return fmt.Errorf("PDBxFileWrite: %w", err)
	}
	defer pdb.Close()
	return PDBxWrite(pdb, coords, mol, bfact, true, strings.Replace(name, ".pdb", "", -1))
}

func PDBxWrite(out io.Writer, coords []*v3.Matrix, mol Atomer, bfact [][]float64, save bool, name ...string) error {
	n := "gochem"
	if len(name) > 0 && name[0] != "" {
		n = name[0]
	}
	out.Write([]byte(fmt.Sprintf("data_%s\n#\n", n)))
	for i, v := range coords {
		cr, _ := v.Dims()
		if cr != mol.Len() {
			return fmt.Errorf("pdbxWrite: Reference (%d) and Coords (%d) don't have the same number of atoms", mol.Len(), cr)
		}
		model := i + 1
		if model == 1 || save {
			out.Write([]byte("loop_\n"))
		}
		if model == 1 {
			out.Write([]byte("_atom_site.type_symbol\n"))
			out.Write([]byte("_atom_site.auth_atom_id\n"))
			out.Write([]byte("_atom_site.auth_comp_id\n"))
			out.Write([]byte("_atom_site.label_alt_id\n"))
			out.Write([]byte("_atom_site.auth_asym_id\n"))
			out.Write([]byte("_atom_site.id\n"))
			out.Write([]byte("_atom_site.auth_seq_id\n"))
			out.Write([]byte("_atom_site.occupancy\n"))
			out.Write([]byte("_atom_site.pdbx_formal_charge\n"))
			out.Write([]byte("_atom_site.group_PDB\n"))
		}
		if model == 1 || save {
			out.Write([]byte("_atom_site.pdbx_PDB_model_num\n"))
			out.Write([]byte("_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"))
			if len(bfact) > i && len(bfact[i]) > mol.Len() {
				out.Write([]byte("_atom_site.B_iso_or_equiv\n"))
			}
		}

		for j := 0; j < mol.Len(); j++ {
			line := ""
			if model == 1 || !save {
				a := mol.Atom(j)
				het := "ATOM"
				if a.Het {
					het = "HETATM"
				}
				c16 := checkChar16(a.Char16)
				line = fmt.Sprintf("%s %s %s %s %s %d %d %4.2f %3.1f %s ", a.Symbol, a.Name, a.MolName, string(c16), a.Chain, a.ID, a.MolID, a.Occupancy, a.Charge, het)
			}
			line += fmt.Sprintf("%d %5.3f %5.3f %5.3f ", model, v.At(j, 0), v.At(j, 1), v.At(j, 2))
			if len(bfact) > i && len(bfact[i]) > mol.Len() {
				line += fmt.Sprintf("%5.3f\n", bfact[i][j])
			} else {
				line += fmt.Sprintf("\n")

			}
			out.Write([]byte(line))
		}
		if save {
			out.Write([]byte("#\n"))
		}
	}
	out.Write([]byte("#\n\n"))
	return nil

}

func checkChar16(c byte) byte {
	valid := "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	for i := range valid {
		if c == valid[i] {
			return c
		}
	}
	return '?'
}

var ma map[string]int = map[string]int{
	"_atom_site.group_pdb":          -1,
	"_atom_site.id":                 -1,
	"_atom_site.type_symbol":        -1,
	"_atom_site.label_atom_id":      -1,
	"_atom_site.label_alt_id":       -1,
	"_atom_site.label_comp_id":      -1,
	"_atom_site.label_asym_id":      -1,
	"_atom_site.label_entity_id":    -1,
	"_atom_site.label_seq_id":       -1,
	"_atom_site.pdbx_pdb_ins_code":  -1,
	"_atom_site.cartn_x":            -1,
	"_atom_site.cartn_y":            -1,
	"_atom_site.cartn_z":            -1,
	"_atom_site.occupancy":          -1,
	"_atom_site.b_iso_or_equiv":     -1,
	"_atom_site.pdbx_formal_charge": -1,
	"_atom_site.auth_seq_id":        -1,
	"_atom_site.auth_comp_id":       -1,
	"_atom_site.auth_asym_id":       -1,
	"_atom_site.auth_atom_id":       -1,
	"_atom_site.pdbx_pdb_model_num": -1,
}
