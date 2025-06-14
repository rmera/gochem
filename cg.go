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

// BackboneCGize takes a coord and a mol for a protein, and returns a new set of coordinates
// each of which cooresponds to the center of mass of the backbone of the corresponding residue
// in the original molecule. If top is true, it also returns a topology where each atom corrsponds
// to the center of mass, with the name "BB", and the correct residue name and ID. Otherwise it returns
// an empty topology.
// In other words, it takes a protein and returns a CG model for its backbone, in the currect conformation.
// If centroid is given and true, the geometric center is used instead of the center of mass.
func BackboneCGize(coord *v3.Matrix, mol Atomer, top bool, centroid ...bool) (*v3.Matrix, *Topology, error) {
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
		var mass []float64
		if len(centroid) > 0 && centroid[0] {
			mass = make([]float64, 0, bbtop.Len())
			for i := 0; i < bbtop.Len(); i++ {
				mass = append(mass, 1)
			}
		} else {
			var err error
			mass, err = bbtop.Masses()
			if err != nil {
				return nil, nil, errorDec(err, i)
			}
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
