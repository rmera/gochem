package gro

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
// from the last atom of F. AType and LJ are merged by preserving all
// elements in the receiver plus the non-repeated elements of F2, if keepours is true
// or preserving all elements in F2 plus the non-repeated elements in the receiver,
// otherwise. Note that 2 LJ elements could be repeated (i.e. refer to the same atoms)
// but have different values, so the choice between keepours true or false is important.
func (F *FF) Merge(F2 *FF, keepours bool) {
	disp := F.Mol.Len()
	sum := func(i int) int {
		return i + disp
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
// it only returns error if an invalid term is given, so, if that is
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
		if !slices.Contains(todel, v.ID) {
			vsites = append(vsites, v)
		}
	}
	F.Mol.Atoms = nat
	F.VSites = vsites
	return F
}
