package gro

import (
	"fmt"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
)

// Returns a string without gromacs comments (sequences starting with ';'),
// trailing and leading spaces, tabs and newlines
func cleanString(s string) string {
	f := strings.Split(s, ";")[0]
	return strings.Trim(f, "\n\t ")

}

// Adds the data in the gromacs-topology atom-section string to the atom with index index[0] in the
// molecule in ff. IF index is not given, the function will search the molecule to add the data to the
// atom that matches the ID on the topology string.
func (F *FF) GroTopAtom2Atom(s string, index ...int) (err error) {
	s = cleanString(s)
	var at *chem.Atom
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	ix := -1
	if len(index) > 0 {
		ix = index[0]
	}
	l := fi(s)
	ID, err := strconv.Atoi(l[0])
	if ix > 0 && F.Mol.Atom(ix).ID == ID {
		at = F.Mol.Atom(ix)
	} else {
		for i := 0; i < F.Mol.Len(); i++ {
			a := F.Mol.Atom(i)
			if a.ID == ID {
				at = a
			}
		}
	}
	if at == nil {
		err = fmt.Errorf("Couldn't find atom with ID %d", ID)
		return
	}
	qerr(err)
	at.Symbol = l[1]
	at.MolID, err = strconv.Atoi(l[2])
	qerr(err)
	at.MolName = l[3]
	at.Name = l[4]
	at.Charge, err = strconv.ParseFloat(l[6], 64)
	qerr(err)
	return nil
}

// Writes the ith (0-based) atom in the FF molecule to a Gromacs sopology line
func (F *FF) Atom2GroTop(i int) string {
	A := F.Mol.Atom(i)
	return fmt.Sprintf("    %5d  %4s %5d  %4s  %4s %5d  %6.4f \n", A.ID, A.Symbol, A.MolID, A.MolName, A.Name, A.ID, A.Charge)
}

// Returns a term containing the information in the GromacsTop-formatted string s,
// given that the string is part of the header header.
func TermFromGroTop(s, header string) (T *Term, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	T = new(Term)
	l := strings.Fields(cleanString(s))
	switch header {
	case "dihedrals":
		if len(l) == 11 {
			// Ryckaert-Bellemans entonces
			T.IDs, err = parseints(l[:4]...)
			qerr(err)
			T.Functype = 3
			T.RB, err = parsefloats(l[5:]...)
			qerr(err)
			if len(T.RB) != 6 {
				err = fmt.Errorf("R-B term detected but read %d parameters instead of the 6 expected", len(T.RB))
				T = nil
				return
			}
			return

		}
	case "constraints":
		//constraints
		T.IDs, err = parseints(l[:2]...)
		T.Constraint = true
		T.Eq, err = strconv.ParseFloat(l[2], 64) //could be 3 if there is a function type but I don't think there is. Check
		qerr(err)
		return

	}
	//Gotta add vsite support

	//now the 'regular' cases

	return
}

func (T *Term) writeAtoms() string {
	add := T.OneBased
	r := make([]string, 0, len(T.IDs))
	for _, v := range T.IDs {
		r = append(r, fmt.Sprintf("%4d", v+add))
	}
	return strings.Join(r, " ")

}

// Writes the term to a string in Gromacs top format.
func (T *Term) ToGroTop() string {
	ret := make([]string, 0, 4)
	ret = append(ret, T.writeAtoms())
	if T.Vsite {
		return "" //placeholder
	}
	if !T.Constraint {
		ret = append(ret, fmt.Sprintf("%1d", T.Functype))
	}
	ret = append(ret, fmt.Sprintf("%6.4f", T.Eq))
	if !T.Constraint && len(T.RB) == 0 {
		ret = append(ret, fmt.Sprintf("%6.4f", T.K))
	}
	if len(T.RB) > 0 && len(T.RB) < 6 {
		panic(fmt.Sprintf("grotop/Term.String: R-B potential must have 6 parameters, got %d", len(T.RB)))
	}
	for _, v := range T.RB {
		ret = append(ret, fmt.Sprintf("%6.4f", v))
	}
	return strings.Join(ret, " ")
}

// Returns a *VSite witht he information in
// the string s containing a virtual_sitesn line
func VSitesNFromGroTop(s string) (vsit *VSite, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	vsit = new(VSite)
	f := strings.Fields(cleanString(s))
	id, err := strconv.Atoi(f[0])
	qerr(err)
	typ, err := strconv.Atoi(f[1])
	qerr(err)
	if id < 0 || typ < 0 {
		return nil, fmt.Errorf("ill-formatted vsiten string: %s", s)
	}
	vsit.ID = id
	vsit.FuncType = typ
	vsit.N = 0 //marks a vsitesn
	var ids []int
	if typ != 3 {
		ids = make([]int, 0, len(f[2:]))
		for i, v := range f[2:] {
			if i >= typ {
				break //we don't try to read too many atoms
			}
			num, err := strconv.Atoi(v)
			if err != nil && num < 0 {
				return nil, fmt.Errorf("ill-formatted vsiten string: %s Error: %w", s, err)
			}
			ids = append(ids, num)
		}
	} else {
		terms := len(f[2:])
		if terms%2 != 0 {
			return nil, fmt.Errorf("ill-formatted vsiten string for site with weights: %s Error: %w", s, err)
		}
		ws := make([]float64, 0, len(f[2:])/2)
		for i := 0; i < (terms - 1); i++ {
			num, err := strconv.Atoi(f[2+i])
			if err != nil && num < 0 {
				return nil, fmt.Errorf("ill-formatted vsiten string: %s Error: %w", s, err)
			}
			weight, err := strconv.ParseFloat(f[3+i], 64)
			if err != nil && num < 0 {
				return nil, fmt.Errorf("ill-formatted vsiten string: %s Error: %w", s, err)
			}
			ids = append(ids, num)
			ws = append(ws, weight)
		}
		vsit.Factors = ws
	}
	vsit.Atoms = ids
	return vsit, err
}

// Returns a Gromacstop-formatted virtual_siteX string with the information in the receiver
// 'X' depends on the information in the receiver (the field 'N').
func (V *VSite) ToGro() (string, error) {
	ret := make([]string, 0, 4)
	ret = append(ret, sf("%5d", V.ID))
	ret = append(ret, sf("%2d", V.FuncType))
	if V.N == 0 { //vsites_n
		if V.FuncType == 3 && len(V.Factors) != len(V.Atoms) {
			return "", fmt.Errorf("virtual_site n with weights detected, but weights don't mach atoms")
		}
		for i, v := range V.Atoms {
			if V.FuncType == 3 {
				ret = append(ret, sf("%5d %5.3f", v, V.Factors[i])) //each couple is an atom and a weight
			} else {
				ret = append(ret, sf("%5d", v)) //each couple is an atom and a weight
			}
		}
	} else {
		for _, v := range V.Atoms {
			ret = append(ret, sf("%5d", v)) //each couple is an atom and a weight
		}
		//NOTE:Right now I'm assuming factors are just in Gromacs units, which
		//is not consistent with 'regular' goChem units.
		//Support for vsites other than 'n' is not quite there yet.
		for _, v := range V.Factors {
			ret = append(ret, sf("%5.3f", v)) //each couple is an atom and a weight
		}

	}
	ret = append(ret, "\n")
	return strings.Join(ret, " "), nil //I'm not sure about the separator. Just "" could be enough.

}
