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
)

//constants from DOI:10.1186/1758-2946-3-33
const (
	tooclose = 0.63
	bondtol  = 0.045
)

type Bond struct {
	Index int
	At1   *Atom
	At2   *Atom
	Dist  float64
	Order float64 //Order 0 means undetermined
}

//Cross returs the atom bonded to the origin atom
//by the received bond.
func (B *Bond) Cross(origin *Atom) *Atom {
	if origin.Index == B.At1.Index {
		return B.At2
	}
	if origin.Index == B.At2.Index {
		return B.At1
	}
	panic("Trying to cross a bond: The origin atom given is not present in the bond!") //I think this got to be a programming error, so a panic is warranted.

}

//Removes the receiver bond from the the Bond slices in the corresponding atoms
func (b *Bond) Remove() error {
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

//ShortestOrLongestPath determines the shortest or longest bonded path between at and the atom with
//Index targetIndex. It returns a slice containing all the atoms in between at and the target (including
//the Index of at) or nil, if there is no valid path between the atoms. "path" is the the path
//already "walked", so, in a new search, it should not be given (although a nil value, or an empty slice
//will also work. All atoms in the molecule need to have the "Index" field filled.
//Even the targetIndex is the same as that of the current atom, and the path is not given, nil, or
//of len 0, the function will not return. This means that the function can be used to search for
//a cyclic path back to the initial atom.
func ShortestOrLongestPath(at *Atom, targetIndex int, shortest bool, path ...[]int) []int {
	//	fmt.Println("yeeeey", at.Index, targetIndex, at.Bonds) /////////////

	if len(path) > 0 && len(path[0]) > 1 && at.Index == path[0][len(path[0])-2] {
		return nil //means we took the same bond back to the previous node, not a valid path.
	}
	if len(path) == 0 {
		path = append(path, []int{at.Index})
	} else if path[0] == nil {
		path[0] = []int{at.Index}

	} else {
		path[0] = append(path[0], at.Index)
	}
	fmt.Printf("last visited is %v\n", path[0])
	if at.Index == targetIndex && len(path[0]) > 1 {
		return path[0] //We arrived! Note that if the starting node is the same as the target, we will
		//only settle for a "cyclic" path that goes through at least another atom (really, at least 2 more atoms).
		// We will not immediately return success on the first node. This is enforced by the len(path[0]>1 condition.
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
