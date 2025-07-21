/*
 * cg.go, part of gochem
 *
 * Copyright 2020 Raul Mera A. (raulpuntomeraatusachpuntocl)
 *
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2.1 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
*/

package chem

import (
	"fmt"
	"log"

	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

func errorDec(err error, frame int) error {
	if err, ok := err.(Error); ok {
		err.Decorate(fmt.Sprintf("BackboneCGize: Error processing %d th residue", frame))
	}
	return err

}

// returns a slice with the masses and an error or nil. If usecom is not given,
// a slice of 1s and nil will be returned (i.e. a slice that can be
// used to calculate the COG).
// if the usecom IS given, but the masses are not found, a slice of 1s
// and an error will be returned.
func masses(top *Topology, usecom ...bool) ([]float64, error) {
	var err error
	var mass []float64
	if len(usecom) > 0 && usecom[0] {
		mass, err = top.Masses()
		if err == nil {
			return mass, nil
		}
	}
	mass = make([]float64, 0, top.Len())
	for i := 0; i < top.Len(); i++ {
		mass = append(mass, 1)
	}
	return mass, err
}

// BackboneCGize takes a coord and a mol for a protein, and returns a new set of coordinates
// each of which cooresponds to the center of mass of the backbone of the corresponding residue
// in the original molecule. If top is true, it also returns a topology where each atom corrsponds
// to the center of mass, with the name "BB", and the correct residue name and ID. Otherwise it returns
// an empty topology.
// In other words, it takes a protein and returns a CG model for its backbone, in the currect conformation.
// If centroid is given and true, the geometric center is used instead of the center of mass.
func BackboneCGize(coord *v3.Matrix, mol Atomer, top bool, com ...bool) (*v3.Matrix, *Topology, error) {
	topol := NewTopology(0, 1)
	residues := countResidues(mol) //sorry
	res4vec := 0
	retcoord := v3.Zeros(residues)
	bb := v3.Zeros(4) //atoms in the backbone
	bbtop := NewTopology(0, 1)
	resid := 0                       //the first residue is 1, not 0!
	for i := 0; i < mol.Len(); i++ { //Note that MolIDs are the PDB residue ID, so they are counted from 1
		atom := mol.Atom(i)
		if resid == atom.MolID {
			continue
		}
		resid = atom.MolID
		bbin := molecules2BBAtoms(mol, []int{resid}, []string{atom.Chain})
		if len(bbin) == 0 {
			continue //not a protein residue, perhaps a ligand
		}
		//	fmt.Println(len(bbin), bbin) ////////////////
		if len(bbin) != bb.NVecs() {
			bb = v3.Zeros(len(bbin))
			if len(bbin) != 4 {
				log.Printf("Abnormal backbone detected in residue %s %d: %d atoms in backbone. Will obtain BB bead with the position of those atoms", atom.MolName, atom.MolID, len(bbin))
			}
		}
		bb.SomeVecs(coord, bbin)
		bbtop.SomeAtoms(mol, bbin)
		if len(com) != 0 {
			com = append(com, true)
		}
		var mass, err = masses(bbtop, com...)
		if err != nil {
			return nil, nil, errorDec(err, i)
		}
		com, err := CenterOfMass(bb, mat.NewDense(len(mass), 1, mass))
		if err != nil {
			return nil, nil, errorDec(err, i)
		}
		retcoord.SetVecs(com, []int{res4vec})
		res4vec++
		if top {
			at := new(Atom)
			at.Copy(atom)
			at.ID = resid
			at.Name = "BB"
			topol.AppendAtom(at)
		}

	}
	if top && topol.Len() != retcoord.NVecs() {
		rc2 := v3.Zeros(topol.Len())
		ri := make([]int, topol.Len())
		for i := 0; i < topol.Len(); i++ {
			ri[i] = i
		}
		rc2.SomeVecs(retcoord, ri)
		retcoord = rc2
	}
	return retcoord, topol, nil
}

func countResidues(mol Atomer) int {
	resid := 0
	ret := 0
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if resid != at.MolID {
			ret++
			resid = at.MolID
		}
	}
	return ret

}

// molecules2BBAtoms returns the indexes for the BB atoms in a given set for residues
// of course, it only works for aminoacidic residues.
func molecules2BBAtoms(mol Atomer, residues []int, chains []string) []int {
	atlist := make([]int, 0, len(residues)*3)
	bb := []string{"N", "CA", "C", "O", "OC1", "OC2"} //OC1 and OC2 are for C-terminal residues which have 2 O.
	for key := 0; key < mol.Len(); key++ {
		at := mol.Atom(key)
		if isInInt(residues, at.MolID) && (isInString(chains, at.Chain) || len(chains) == 0) {
			if isInString(bb, at.Name) {
				atlist = append(atlist, key)
			}
		}
	}
	return atlist
}

// Returns COM or COG (default) coordinates for the beads constructed each from the atoms
// with indexes (in mol) given in the respective slice of beads, and with coordinates coords.
func At2CGCoords(coords *v3.Matrix, mol *Topology, beads [][]int, usecom ...bool) (*v3.Matrix, error) {
	ret := v3.Zeros(len(beads))
	for i, v := range beads {
		tmp := v3.Zeros(len(v))
		tmpt := NewTopology(0, 1)
		tmp.SomeVecs(coords, v)
		tmpt.SomeAtoms(mol, v)
		mass, err := masses(tmpt, usecom...)
		if err != nil {
			return nil, fmt.Errorf("Error assigning masses: %w", err)
		}
		com, err := CenterOfMass(tmp, mat.NewDense(len(mass), 1, mass))
		if err != nil {
			return nil, fmt.Errorf("Couldn't obtain the COM or centroid: %w", err)
		}
		ret.SetVecs(com, []int{i})
	}
	return ret, nil
}
