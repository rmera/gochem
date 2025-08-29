package top

import (
	"fmt"
	"slices"

	chem "github.com/rmera/gochem"
)

func applyToTerms(f func(int) int, T []*Term) {
	for i, v := range T {
		for j, w := range v.IDs {
			T[i].IDs[j] = f(w)
		}
	}
}

// Modifies the IDs of all atoms in the atoms, terms
// and vsite definitions according to f.
func (F *FF) ModIDs(f func(int) int) {
	for i := 0; i < F.Mol.Len(); i++ {
		at := F.Mol.Atom(i)
		at.SetIndex(f(at.Index()))
		at.ID = f(at.ID)
	}
	applyToTerms(f, F.Bonds)
	applyToTerms(f, F.Impropers)
	applyToTerms(f, F.Constraints)
	applyToTerms(f, F.Angles)

	applyToTerms(f, F.Dihedrals)

	for _, v := range F.Exclusions {
		for j, w := range v {
			v[j] = f(w)
		}
	}
	for _, v := range F.VSites {
		v.ID = f(v.ID)
		for j, w := range v.Atoms {
			v.Atoms[j] = f(w)
		}
	}

}

func mergeTerms[T ~[]E, E any](T1, T2 T) T {
	ret := make(T, 0, len(T2)+len(T1))
	for _, v := range T1 {
		ret = append(ret, v)
	}
	for _, v := range T2 {
		ret = append(ret, v)
	}
	return ret

}

// Merges F2 into F, and _after_ F. F2 is modified, as all atom indexes are set to start
// from the last atom of F, unless moveindexes is given and false. AType and LJ are
// merged by preserving all elements in the receiver plus the non-repeated elements of
// F2, if keepours is true or preserving all elements in F2 plus the non-repeated
// elements in the receiver, otherwise. Note that 2 LJ elements could be repeated
// (i.e. refer to the same atoms) but have different values, so the choice between
// keepours true or false is important.
func (F *FF) Merge(F2 *FF, keepours bool, moveindexes ...bool) {
	disp := F.Mol.Len()
	sum := func(i int) int {
		return i + disp
	}
	if len(moveindexes) > 0 || !moveindexes[0] {
		sum = func(i int) int {
			return i
		}

	}
	if keepours {
		F2.ATypes = deleteRepeated(F.ATypes, F2.ATypes)
		F2.LJ = deleteRepeated(F.LJ, F2.LJ)
	} else {
		F.ATypes = deleteRepeated(F2.ATypes, F.ATypes)
		F.LJ = deleteRepeated(F2.LJ, F.LJ)

	}
	F2.ModIDs(sum)
	for i := 0; i < F2.Mol.Len(); i++ {
		F.Mol.AppendAtom(F2.Mol.Atom(i))
	}
	F.Bonds = mergeTerms(F.Bonds, F2.Bonds)
	F.Constraints = mergeTerms(F.Constraints, F2.Constraints)
	F.Angles = mergeTerms(F.Angles, F2.Angles)
	F.Impropers = mergeTerms(F.Impropers, F2.Impropers)
	F.Dihedrals = mergeTerms(F.Dihedrals, F2.Dihedrals)
	F.VSites = mergeTerms(F.VSites, F2.VSites)
	F.Exclusions = mergeTerms(F.Exclusions, F2.Exclusions)
}

type eq interface {
	equal(any) bool
}

func containsEqual(w []eq, e eq) bool {
	for _, v := range w {
		if e.equal(v) {
			return true
		}
	}
	return false

}

// Returns a copy of D whitout the elements that are present in N
func deleteRepeated[E eq, T []E](N, D T) T {
	ret := make([]E, 0, len(N))
	todelete := make([]int, 0, 100)
	for i, r := range D {
		for _, k := range N {
			if r.equal(k) {
				todelete = append(todelete, i)
				break
			}
		}
	}
	for i, v := range D {
		if !slices.Contains(todelete, i) {
			ret = append(ret, v)
		}
	}
	return ret
}

func containsSame[T ~[]E, E comparable](C1, C2 T) bool {
	if len(C1) != len(C2) {
		return false
	}
	for i, v := range C1 {
		if v == C2[i] {
			continue
		}
		r := false
		for _, w := range C2 {
			if v == w {
				r = true
			}
		}
		if !r {
			return r
		}

	}
	return true

}

// Find the Terms with IDs that match the atom IDs in IDs. The match comparison
// depends on the value of 'order'. Acceptable values are:
// 's' (strict, only a term with the values in the same order as
// given is considered).
// 'r' (reverse, both the given order and its reverse are accepted)
// 'a' (any, any order is accepted)
// it only returns error if an invalid order is given, so, if that is
// hardcoded, you may safely omit the error check.
func FindTerm(T []*Term, IDs []int, order byte) (int, error) {
	j := -1
	rev := make([]int, 0, len(IDs))
	copy(rev, IDs)
	slices.Reverse(rev)
	for i, v := range T {
		switch order {
		case 's':
			if slices.Equal(v.IDs, IDs) {
				j = i
				break
			}
		case 'r':
			if slices.Equal(v.IDs, IDs) || slices.Equal(v.IDs, rev) {
				j = i
				break
			}
		case 'a':
			if containsSame(v.IDs, IDs) {
				j = i
			}
		default:
			return -1, fmt.Errorf("Invalid value for order: %v", order)

		}
	}
	return j, nil
}

// Deletes all atoms (and, if applicable, their virtual site definitions)
// with IDs present in todel. It modifies F in place, and also returns it.
// It does _not_ delete terms associated with the atoms.
func DeleteAtomsAndVSites(F *FF, todel []int) *FF {
	//First we generate a new list of atoms without the ones to delete
	nat := make([]*chem.Atom, 0, F.Len()-len(todel))
	for i := 0; i < F.Len(); i++ {
		if !slices.Contains(todel, i+1) {
			nat = append(nat, F.Mol.Atom(i))
		}
	}
	//Now the same for the VSites
	vsites := make([]*VSite, 0, len(F.VSites)) //we don't know how many vsites we'll delete
	for _, v := range F.VSites {
		//	fmt.Println(v.ID, todel) ///////////////////////
		if !slices.Contains(todel, v.ID) {
			vsites = append(vsites, v)
		}
	}
	F.Mol.Atoms = nat
	F.VSites = vsites
	return F
}

// Returns FF like F but removing all atoms, bonded terms and exclusions that _only_
// contain terms in todel or, if removecontaining[0] is given and true, all
// bonded terms and exclusions that contains at least one of the terms in todel.
func DeleteTermsForAtomsAndVSites(F *FF, todel []int, removeallcontaining ...bool) *FF {
	var dodel func([]int) bool
	//here we decide what to delete
	dodel = func(ids []int) bool {
		//this one only deletes if ids ONLY contains atoms that
		//are in todel
		if len(ids) > len(todel) {
			return false
		}
		for _, v := range ids {
			if !slices.Contains(todel, v) {
				return false
			}
		}
		return true
	}
	if len(removeallcontaining) > 0 && removeallcontaining[0] {
		dodel = func(ids []int) bool {
			//this one only deletes if ids has _any_ atom in todel
			for _, v := range ids {
				if slices.Contains(todel, v) {
					return true
				}
			}
			return false
		}

	}
	//now a function to delete terms, since we'll do that a few times
	var DelTerms func([]*Term) []*Term = func(ts []*Term) []*Term {
		ret := make([]*Term, 0, len(ts))
		for _, v := range ts {
			if dodel(v.IDs) {
				continue
			}
			t2 := new(Term)
			t2.Copy(v)
			ret = append(ret, t2)

		}
		return ret
	}

	//Well, that was a lot, but now just sit back
	//and enjoy.

	//Copy the basic data. It's all deep copy as we dont' want
	//to alter the original.
	mol := chem.NewTopology(0, 1)
	mol.CopyAtoms(F.Mol)
	F2 := NewFF(mol, F.SigmaEpsilon)
	F2.currentHeader = F.currentHeader
	F2.VSites = F.VSites
	F2 = DeleteAtomsAndVSites(F2, todel)

	//Now we just apply the functions we worked hard for
	//above.
	F2.Bonds = DelTerms(F.Bonds)
	F2.Angles = DelTerms(F.Angles)
	F2.Constraints = DelTerms(F.Constraints)
	F2.Impropers = DelTerms(F.Impropers)
	F2.Dihedrals = DelTerms(F.Dihedrals)

	//And the exclusions. We still get to recycle the dodel function.
	F2.Exclusions = make([][]int, 0, len(F.Exclusions))
	for _, v := range F.Exclusions {
		if dodel(v) {
			continue
		}
		ex := make([]int, len(v))
		copy(ex, v)
		F2.Exclusions = append(F2.Exclusions, ex)
	}

	return F2
}

func termHole(T []*Term, after, size int) []*Term {
	for i, v := range T {
		for j, w := range v.IDs {
			if w > after {
				T[i].IDs[j] += size
			}
		}
	}
	return T
}

// Returns FF like F but removing all atoms, bonded terms and exclusions that _only_
// contain terms in todel or, if removecontaining[0] is given and true, all
// bonded terms and exclusions that contains at least one of the terms in todel.
func Shift(F *FF, shift int, shiftAtoms bool) {
	//here we decide what to delete
	//now a function to delete terms, since we'll do that a few times
	var ShiftTerms func([]*Term) = func(ts []*Term) {
		for _, v := range ts {
			for j, _ := range v.IDs {
				v.IDs[j] += shift
			}

		}
	}
	ShiftTerms(F.Bonds)
	ShiftTerms(F.Angles)
	ShiftTerms(F.Constraints)
	ShiftTerms(F.Impropers)
	ShiftTerms(F.Dihedrals)
	for i, v := range F.Exclusions {
		for j, _ := range v {
			F.Exclusions[i][j] += shift
		}
	}

	if shiftAtoms {
		for i := 0; i < F.Len(); i++ {
			at := F.Mol.Atom(i)
			at.ID += shift
			at.SetIndex(at.Index() + shift)
		}
	}
	//Vsites get shifted even if atoms don't.
	for i, _ := range F.VSites {
		F.VSites[i].ID += shift

	}
}

// Changes all ids in the given terms according to the given map
// modifies the terms and returns them.
func switchterms(m map[int]int, te []*Term) []*Term {
	for i, v := range te {
		for j, w := range v.IDs {
			te[i].IDs[j] = m[w]
		}
	}
	return te
}

func SwitchMap(oripos, newpos []int, fulllen int) map[int]int {
	m := make(map[int]int)
	for i := 0; i < fulllen; i++ {
		m[i] = i
	}
	for i, v := range oripos {
		m[v] = newpos[i]
	}
	return m
}

// Switches atoms from originalposition to newposition
func Switch(originalposition, newposition []int, F *FF) *FF {
	o := originalposition
	n := newposition
	m := SwitchMap(o, n, F.Mol.Len())
	switchterms(m, F.Bonds)
	switchterms(m, F.Angles)
	switchterms(m, F.Constraints)
	switchterms(m, F.Impropers)
	switchterms(m, F.Dihedrals)
	F.Mol = chem.SwitchAtoms(o, n, F.Mol)
	for i, v := range F.Exclusions {
		for j, w := range v {
			F.Exclusions[i][j] = m[w]
		}
	}
	//Vsites get shifted even if atoms don't.
	for i, _ := range F.VSites {
		F.VSites[i].ID += m[F.VSites[i].ID]
	}
	return F
}

// not ready
/*
func MakeHole(after int, size int, F *FF) {
	F.Bonds = termHole(F.Bonds, after, size)
	F.Angles = termHole(F.Angles, after, size)

	for i, v := range F.Exclusions {
		for j, w := range v {
			if w > after {
				F.Exclusions[i][j] += size
			}
		}
	}

	for i, v := range F.VSites {
		if v.ID > after {
			F.VSites[i].ID += size
		}
		for j, w := range v.Atoms {
			if w > after {
				F.Bonds[i].IDs[j] += size
			}
		}
	}

}
*/

/*
type FF struct {
	SigmaEpsilon  bool //are LJ terms using sigma/epsilon, or C6/C12?
	currentHeader string
	Mol           *chem.Topology
	Bonds         []*Term
	Constraints   []*Term
	Angles        []*Term
	Impropers     []*Term
	Dihedrals     []*Term
	VSites        []*VSite
	ATypes        []*AtomType
	LJ            []*LJPair
	Exclusions    [][]int
}
*/
