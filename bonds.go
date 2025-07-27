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
 *
 */

package chem

import (
	"fmt"
	"slices"
	"sort"

	v3 "github.com/rmera/gochem/v3"
)

// constants from DOI:10.1186/1758-2946-3-33
const (
	tooclose = 0.63
	bondtol  = 0.045 //the reference actually says 0.45 A but I strongly suspect that is a typo.
)

// Bond represents a chemical, covalent bond.
type Bond struct {
	Index  int
	At1    *Atom
	At2    *Atom
	Dist   float64
	Energy float64 //One might prefer to use energy insteaf of order.
	//NOTE
	//I'm not sure if I leave just the order and let the user "Abuse" that field for energy
	//or anything else they want, or if I keep this extra field.
	Order float64 //Order 0 means undetermined
}

// Copies B into the receiver
func (b *Bond) Copy(B *Bond) {
	b.Index = B.Index
	b.At1 = B.At1
	b.At2 = B.At2
	b.Dist = B.Dist
	b.Energy = B.Energy //One might prefer to use energy insteaf of order.
	//NOTE
	//I'm not sure if I leave just the order and let the user "Abuse" that field for energy
	//or anything else they want, or if I keep this extra field.
	b.Order = B.Order //Order 0 means undetermined
}

// Cross returns the atom bonded to the origin atom
// bond in the receiver.
func (B *Bond) Cross(origin *Atom) *Atom {
	if origin.index == B.At1.index {
		return B.At2
	}
	if origin.index == B.At2.index {
		return B.At1
	}
	panic("Trying to cross a bond: The origin atom given is not present in the bond!") //I think this got to be a programming error, so a panic is warranted.

}

// Remove removes the receiver bond from the the Bond slices in the corresponding atoms
func (b *Bond) Remove() error {
	lenb1 := len(b.At1.Bonds)
	lenb2 := len(b.At2.Bonds)
	b.At1.Bonds = takefromslice(b.At1.Bonds, b.Index)
	b.At2.Bonds = takefromslice(b.At2.Bonds, b.Index)
	err := new(CError)
	errs := 0
	err.msg = fmt.Sprintf("Failed to remove bond Index:%d", b.Index)
	if len(b.At1.Bonds) == lenb1 {
		err.msg = err.msg + fmt.Sprintf("from atom. Index:%d", b.At1.index)
		err.Decorate("RemoveBond")
		errs++
	}
	if len(b.At2.Bonds) == lenb2 {
		err := new(CError)
		if errs > 0 {
			err.msg = err.msg + " and"
		}
		err.msg = err.msg + fmt.Sprintf("from atom. Index:%d", b.At2.index)
		err.Decorate("RemoveBond")
		errs++
	}
	if errs > 0 {
		return err
	}
	return nil
}

// returns a new *Bond slice with the element id removed
func takefromslice(bonds []*Bond, id int) []*Bond {
	newb := make([]*Bond, len(bonds)-1)
	for _, v := range bonds {
		if v.Index != id {
			newb = append(newb, v)
		}
	}
	return newb
}

// BondedOptions contains options for the BondePaths function
type BondedOptions struct {
	OnlyShortest bool             //Only return the shortest path between the atoms
	F            func(*Bond) bool //Only bonds for which this function returns true are considered
	path         []int            //
}

func DefaultBondedOptions() *BondedOptions {
	return &BondedOptions{OnlyShortest: false, F: func(*Bond) bool { return true }}
}

/******
**** The following is commented out and excluded from the API, as I don't think it is needed.
**** In any case, adding this function, or any other, wouldn't break the API, so it's best
**** to exclude when in doubt. Of course this is in itself an API break, but the bond functionality
**** is still very new, so the break is very unlikely to affect anyone.
//SetAlreadyWalkedPath set the unexported field "path" in BondedOptions to p.
//the path field represents the already-walked path, so you almost never
//want to set it. The exception is definidng your own function that calls
//BondedPath recursively, as I do here. This method was added to allow
//for such use.
func (B *BondedOptions)SetAlreadyWalkedPath(p []int){
    B.path=p
}
******/

/*
// BondedPaths determines the paths between at and the atom with
// Index targetIndex. It returns a slice of slices of int, where each sub-slice contains all the atoms
// in between at and the target (including the index of at)
// If there is no valid path, it returns nil.
// In the options structure BondeOption, OnlyShortest set to true causes only the shortest of the found paths to
// be returned.  All atoms in the molecule need to have the "index" field filled.
// If the targetIndex is the same as that of the current atom, and the path is not given, nil, or
// of len 0, the function will search for a cyclic path back to the initial atom.
// if onlyshortest is true, only the shortest path will be returned (the other elements of the slice will be nil)
// This can be useful if you want to save memory on a very intrincate molecule.
func BondedPaths(at *Atom, targetIndex int, options ...*BondedOptions) [][]int {
	if len(options) == 0 {
		options = []*BondedOptions{{OnlyShortest: false, path: nil}}
	}
	onlyshortest := options[0].OnlyShortest
	path := [][]int{options[0].path}
	//I am not completely sure about this function signature. It is a candidate for API change.
	if len(path) > 0 && len(path[0]) > 1 && path[0][len(path[0])-2] == at.index {
		return nil //We are back to the atom we just had visited, not a valid path. We have to check this before checking if we completed the "quest"
		//or, by just going back via the same bond, it would seem like we are at the finishing line.
	}
	if len(path) == 0 {
		path = append(path, []int{at.index})
	} else if path[0] == nil {
		path[0] = []int{at.index}

	} else {
		path[0] = append(path[0], at.index)
	}

	if at.index == targetIndex && len(path[0]) > 1 {
		return [][]int{path[0]} //We arrived! Note that if the starting node is the same as the target, we will
		//only settle for a "cyclic" path that goes through at least another atom (really, at least 2 more atoms).
		// We will not immediately return success on the first node. This is enforced by the len(path[0]>1 condition.
	}
	//Here we check that we are not back to an atom we previously visited. This checks for loops, and has to be performed
	//after we check if we got to the goal, since the goal could be the same starting atom (if we look for a cyclic path).
	if len(path[0]) > 1 && isInInt(path[0][:len(path[0])-1], at.index) {
		return nil //means we took the same bond back to the previous node, or got trapped in a loop. not a valid path.
	}
	if len(at.Bonds) <= 1 {
		return nil //means that we hit an end of the road. There is only one bond in the atom (i.e. the one leading to the previous node)
	}
	rets := make([][]int, 0, len(at.Bonds))
	for _, v := range at.Bonds {
		path2 := make([]int, len(path[0]))
		copy(path2, path[0])
		rets = append(rets, BondedPaths(v.Cross(at), targetIndex, &BondedOptions{OnlyShortest: onlyshortest, path: path2})...) //scary stuff
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
	if onlyshortest {
		return [][]int{rets2[0]}
	}
	return rets2

}
*/

// BondedPaths determines the paths between at and the atom with
// Index targetIndex. It returns a slice of slices of int, where each sub-slice contains all the atoms
// in between at and the target (including the index of at)
// If there is no valid path, it returns nil.
// In the options structure BondeOption, OnlyShortest set to true causes only the shortest of the found paths to
// be returned.  All atoms in the molecule need to have the "index" field filled.
// If the targetIndex is the same as that of the current atom, and the path is not given, nil, or
// of len 0, the function will search for a cyclic path back to the initial atom.
// if onlyshortest is true, only the shortest path will be returned (the other elements of the slice will be nil)
// This can be useful if you want to save memory on a very intrincate molecule.
// Only bonds for which the function in options returns true are considered for the search (in the default options,
// this function always returns true, so all bonds are consiered).
func BondedPaths(at *Atom, targetIndex int, options ...*BondedOptions) [][]int {
	if len(options) == 0 {
		options = append(options, DefaultBondedOptions())
	}
	onlyshortest := options[0].OnlyShortest
	path := [][]int{options[0].path}
	if len(path) > 0 && len(path[0]) > 1 && path[0][len(path[0])-2] == at.index {
		return nil //We are back to the atom we just had visited, not a valid path. We have to check this before checking if we completed the "quest"
		//or, by just going back via the same bond, it would seem like we are at the finishing line.
		//This really shouldn't happen, but I'm paranoid.
	}
	if len(path) == 0 {
		path = append(path, []int{at.index})
	} else if path[0] == nil {
		path[0] = []int{at.index}

	} else {
		path[0] = append(path[0], at.index)
	}

	if at.index == targetIndex && len(path[0]) > 1 {
		return [][]int{path[0]} //We arrived! Note that if the starting node is the same as the target, we will
		//only settle for a "cyclic" path that goes through at least another atom (really, at least 2 more atoms).
		// We will not immediately return success on the first node. This is enforced by the len(path[0]>1 condition.
	}
	//Here we check that we are not back to an atom we previously visited. This checks for loops, and has to be performed
	//after we check if we got to the goal, since the goal could be the same starting atom (if we look for a cyclic path).
	if len(path[0]) > 1 && isInInt(path[0][:len(path[0])-1], at.index) {
		return nil //means we took the same bond back to the previous node, or got trapped in a loop. not a valid path.
	}
	if len(at.Bonds) <= 1 {
		return nil //means that we hit an end of the road. There is only one bond in the atom (i.e. the one leading to the previous node)
	}
	rets := make([][]int, 0, len(at.Bonds))
	for _, v := range at.Bonds {
		if !options[0].F(v) {
			continue
		}
		path2 := make([]int, len(path[0]))
		copy(path2, path[0])
		rets = append(rets, BondedPaths(v.Cross(at), targetIndex, &BondedOptions{OnlyShortest: onlyshortest, path: path2, F: options[0].F})...) //scary stuff
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
	if onlyshortest {
		return [][]int{rets2[0]}
	}
	return rets2

}

// Ring represents a molecular cycle.
type Ring struct {
	Atoms     []int
	planarity float64
}

// Contains returns true or false depending on whether
// the atom with the given index is part of the ring
func (R *Ring) Contains(index int) bool {
	return isInInt(R.Atoms, index)
}

func (R *Ring) IsIn(index int) bool {
	return R.Contains(index)
}

// Size returns the number of atoms in the ring
func (R *Ring) Size() int {
	return len(R.Atoms)
}

// Planarity returns the planarity percentage of the receiver ring
// coords is the set of coordinates for the _entire_ molecule of which
// the ring is part. Planarity does not check that coords indeed corresponds
// to the correct molecule, so, doing so is the user's responsibility.
func (R *Ring) Planarity(coord *v3.Matrix, temp ...*v3.Matrix) float64 {
	if R.planarity != 0 {
		return R.planarity
	}
	var c *v3.Matrix
	if len(temp) > 0 && temp[0] != nil && temp[0].Len() == len(R.Atoms) {
		c = temp[0]
	} else {
		c = v3.Zeros(len(R.Atoms))
	}
	c.SomeVecs(coord, R.Atoms)
	_, plan, err := EasyShape(c, 0.01)
	if err != nil {
		R.planarity = -1
		return -1
	}
	R.planarity = plan
	return plan

}

// AddHs Adds to the ring the H atoms bonded to its members
// mol is the _entire_ molecule of which the receiver ring is part.
// It will panic if an atom of the ring is not
// found in mol, or is not bonded to anything
func (R *Ring) AddHs(mol Atomer) {
	newind := make([]int, 0, 6)
	for _, v := range R.Atoms {
		at := mol.Atom(v) //this can panic if you give the wrong mol object
		for _, w := range at.Bonds {
			at2 := w.Cross(at)
			if at2.Symbol == "H" {
				newind = append(newind, at2.index)
			}

		}
	}
	R.Atoms = append(R.Atoms, newind...)
}

// API CHANGE
// InWhichRing returns the indexes of all rings top which at belongs
// or nil if it doesn't belong to any ring.
func InWhichRing(at *Atom, rings []*Ring) []int {
	if len(rings) == 0 {
		return nil
	}
	var ret []int
	for i, v := range rings {
		if v.IsIn(at.index) {
			ret = append(ret, i)
		}
	}
	return ret

}

type RingOptions struct {
	AddHs        bool
	MinPlanarity float64
	FusedRings   bool
	MaxRing      int //maximum ring size
	Coords       *v3.Matrix
	F            func(*Bond) bool
}

func DefaultRingOptions() *RingOptions {
	f := func(*Bond) bool { return true }
	return &RingOptions{MinPlanarity: -1, Coords: nil, AddHs: false, MaxRing: 6, F: f}

}

// Identifies and returns all rings in mol, by
// searching for cyclic paths. if at least 1 element is given for addHs
// and the first element of those give is true, the Hs bound to the atoms in the rings
// will also be added to the ring.
func FindRings(mol Atomer, Opts ...*RingOptions) []*Ring {
	var o *RingOptions
	if len(Opts) > 0 {
		o = Opts[0]
	} else {
		o = DefaultRingOptions()
	}
	L := mol.Len()
	var rings []*Ring
	for i := 0; i < L; i++ {
		at := mol.Atom(i)
		//Here I nullified the requirement that at is not in a previous ring, if fused rings are wanted.
		//This is because, say you have 2 fused rings, A and B, with the 'metaring', the fusing of both, named C.
		//you detect first A and its 'metaring' C, but then you don't detect B as all its atoms are part of C,
		//and thys gets excluded byt the (InWhichRing(at, rings) == nil) condition
		if o.FusedRings || (InWhichRing(at, rings) == nil) {
			paths := BondedPaths(at, at.index, &BondedOptions{OnlyShortest: !o.FusedRings, F: o.F})
			if len(paths) == 0 || (len(paths[0][1:]) > o.MaxRing && o.MaxRing > 0) {
				continue
			}

			//NOTE: I added an additional loop here
			//to allow for 'fused rings' i.e. longer paths connecting the same atoms.
			//it should only matter if o.FusedRings is set to true, as, otherwise, len(paths)==1
			for _, p := range paths {
				r := &Ring{Atoms: p[1:]} //NOTE: it was originally paths[0], so, if something breaks
				//look here. I think the previous behavior was buggy, since the first atom of the ring will be repeated
				//as the last atom, giving a wrong number for the ring's length.
				if o.Coords == nil || o.MinPlanarity < 0 || r.Planarity(o.Coords) > o.MinPlanarity {
					if o.AddHs {
						r.AddHs(mol)
					}
					rings = append(rings, r)
				}
			}
		}

	}

	return removeRepeatedRings(rings)
}

// 2 rings are equal if they contain the same atoms, whether in the same
// order, or not.
func ringIsEqual(a, b *Ring) bool {
	if len(a.Atoms) != len(b.Atoms) {
		return false
	}
	hasall := true
	for _, v := range a.Atoms {
		if !slices.Contains(b.Atoms, v) {
			hasall = false
			break
		}
	}
	return hasall
}

func removeRepeatedRings(set []*Ring) []*Ring {
	ret := make([]*Ring, 0, 1)
	for _, v := range set {
		rep := false
		for _, w := range ret {
			if ringIsEqual(v, w) {
				rep = true
				break
			}
		}
		if !rep {
			ret = append(ret, v)
		}
	}
	return ret
}
