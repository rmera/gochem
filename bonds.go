/*
 * atomicdata.go, part of gochem.
 *
 *
 * Copyright 2021 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 * goChem is currently developed at the Universidad de Santiago de Chile
 * (USACH)
 *
 */

package chem

import (
	"fmt"
	"sort"

	v3 "github.com/rmera/gochem/v3"
)

//constants from DOI:10.1186/1758-2946-3-33
const (
	tooclose = 0.63
	bondtol  = 0.45
	toofar   = 3
)

type Bond struct {
	Index int
	At1   *Atom
	At2   *Atom
	Dist  float64
	Order float64 //Order 0 means undetermined
}

//Cross takes
func (B *Bond) Cross(origin *Atom) *Atom {
	if origin.Index == B.At1.Index {
		return B.At2
	}
	if origin.Index == B.At2.Index {
		return B.At1
	}
	panic("Trying to cross a bond: The origin atom given is not present in the bond!") //I think this got to be a programming error, so a panic is warranted.

}

//return a new *Bond slice with the element id removed
func takefromslice(bonds []*Bond, id int) []*Bond {
	newb := make([]*Bond, len(bonds)-1)
	for _, v := range bonds {
		if v.Index != id {
			newb = append(newb, v)
		}
	}
	return newb
}

func RemoveBond(b *Bond, mol Atomer) error {
	lenb1 := len(b.At1.Bonds)
	lenb2 := len(b.At2.Bonds)
	b.At1.Bonds = takefromslice(b.At1.Bonds, b.Index)
	b.At2.Bonds = takefromslice(b.At2.Bonds, b.Index)
	err := new(CError)
	errs := 0
	err.msg = fmt.Sprintf("Failed to remove bond Index:%d", b.Index)
	if len(b.At1.Bonds) == lenb1 {
		err.msg = err.msg + fmt.Sprintf("from atom. Index:%d", b.At1.Index)
		err.Decorate("RemoveBond")
		errs++
	}
	if len(b.At2.Bonds) == lenb2 {
		err := new(CError)
		if errs > 0 {
			err.msg = err.msg + " and"
		}
		err.msg = err.msg + fmt.Sprintf("from atom. Index:%d", b.At2.Index)
		err.Decorate("RemoveBond")
		errs++
	}
	if errs > 0 {
		return err
	}
	return nil
}

//Assigns bonds to a molecule based on a simple distance
//criterium, similar to that described in DOI:10.1186/1758-2946-3-33
func AssignBonds(coord *v3.Matrix, mol AtomIndexesFiller) error {
	// might get slow for
	//large systems. It's really not thought
	//for proteins or macromolecules.
	var t1, t2 *v3.Matrix
	var at1, at2 *Atom
	mol.FillIndexes()
	t3 := v3.Zeros(1)
	bonds := make([]*Bond, 0, 10)
	tot := mol.Len()
	var nextIndex int
	for i := 0; i < tot; i++ {
		t1 = coord.VecView(i)
		at1 = mol.Atom(i)
		cov1 := symbolCovrad[at1.Symbol]
		if cov1 == 0 {
			err := new(CError)
			err.msg = fmt.Sprintf("Couldn't find the covalent radii  for %s %d", at1.Symbol, i)
			err.Decorate("AssignBonds")
			return err
		}
		for j := i + 1; j < tot; j++ {
			t2 = coord.VecView(j)
			at2 = mol.Atom(j)
			cov2 := symbolCovrad[at2.Symbol]
			if cov2 == 0 {
				err := new(CError)
				err.msg = fmt.Sprintf("Couldn't find the covalent radii  for %s %d", at2.Symbol, j)
				err.Decorate("AssignBonds")
				return err
			}

			t3.Sub(t2, t1)
			d := t3.Norm(2)
			if d < cov1+cov2+bondtol && d > tooclose {
				b := &Bond{Index: nextIndex, Dist: d, At1: at1, At2: at2}
				at1.Bonds = append(at1.Bonds, b)
				at2.Bonds = append(at2.Bonds, b)
				bonds = append(bonds, b) //just to easily keep track of them.
				nextIndex++
			}

		}
	}

	//Now we check that no atom has too many bonds.
	for i := 0; i < tot; i++ {
		at := mol.Atom(i)
		max := symbolMaxBonds[at.Symbol]
		if max == 0 { //means there is not a specified number of bonds for this atom.
			continue
		}
		sort.Slice(at.Bonds, func(i, j int) bool { return at.Bonds[i].Dist < at.Bonds[j].Dist })
		//I am hoping this will remove bonds until len(at.Bonds) is not
		//greater than max.
		for i := len(at.Bonds); i > max; i = len(at.Bonds) {
			err := RemoveBond(at.Bonds[len(at.Bonds)-1], mol) //we remove the longest bond
			if err != nil {
				return errDecorate(err, "AssignBonds")
			}
		}

	}

	return nil
}

//ShortestOrLongestPath determines the shortest or longest path between at and the atom with
//Index targetIndex. It returns a slice containing all the atoms in between at and the target (including
//the Index of at) or nil, if there is no valid path between the atoms. "path" is the the path
//already "walked", so, in a new search, it should not be given (although a nil value, or an empty slice
//will also work. All atoms in the molecule need to have the "Index" field filled.
//Even the targetIndex is the same as that of the current atom, and the path is not given, nil, or
//of len 0, the function will not return. This means that the function can be used to search for
//a cyclic path back to the initial atom.
func ShortestOrLongestPath(at *Atom, targetIndex int, shortest bool, path ...[]int) []int {
	if len(path) != 0 && path[0] != nil && len(path[0]) > 0 && at.Index == path[0][len(path)-1] {
		return nil //means we took the same bond back to the previous node, not a valid path.
	}
	if len(path) == 0 {
		path = append(path, []int{at.Index})
	} else if path[0] == nil {
		path[0] = []int{at.Index}

	} else {
		path[0] = append(path[0], at.Index)
	}

	if at.Index == targetIndex && len(path[0]) > 1 {
		return path[0] //We arrived! Note that if the starting node is the same as the target, we will
		//only settle for a "cyclic" ring. We will not immediately return success on the first node.
	}
	if len(at.Bonds) <= 1 {
		return nil //means that we hit an end of the road. There is only one bond in the atom (i.e. the one leading to the previous node)
	}
	rets := make([][]int, 0, len(at.Bonds))
	for _, v := range at.Bonds {
		path2 := make([]int, len(path[0]))
		copy(path2, path[0])
		rets = append(rets, ShortestOrLongestPath(v.Cross(at), targetIndex, shortest, path2)) //scary stuff
	}
	rets2 := make([][]int, 0, len(at.Bonds))
	for _, v := range rets {
		if v != nil {
			rets2 = append(rets2, v)
		}
	}
	if len(rets2) == 0 {
		return nil
	}
	sort.Slice(rets2, func(i, j int) bool { return len(rets2[i]) < len(rets2[j]) })
	if shortest {
		return rets2[0]
	}
	return rets2[len(rets2)-1]

}
