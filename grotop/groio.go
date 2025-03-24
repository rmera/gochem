package gro

import (
	"fmt"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
)

func GroTopAtom2Atom(s string, mol chem.Atomer, index ...int) (err error) {
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
	if ix > 0 && mol.Atom(ix).ID == ID {
		at = mol.Atom(ix)
	} else {
		for i := 0; i < mol.Len(); i++ {
			a := mol.Atom(i)
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

func Atom2GroTop(mol chem.Atomer, i int) string {
	A := mol.Atom(i)
	return fmt.Sprintf("    %5d  %4s %5d  %4s  %4s %5d  %6.4f \n", A.ID, A.Symbol, A.MolID, A.MolName, A.Name, A.ID, A.Charge)
}

func TermFromGroTop(s, header string) (T *Term, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	T = new(Term)
	l := strings.Fields(delcomment(s))
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
