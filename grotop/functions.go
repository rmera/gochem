package gro

import (
	"errors"
	"fmt"
	"io"
	"slices"
	"strconv"
	"strings"
)

var fi func(string) []string = strings.Fields
var sf func(string, ...any) string = fmt.Sprintf

// Returns the same function that will delete lines that have been already
// seen in the file. This will use a fair bit of memory as all lines have to be
// kept in memory for comparison with the next ones.
func DelRepeatedFunctions(T TermSelect) {
	lines := make([]string, 0, 200)
	f := func(s string) (string, error) {
		if slices.Contains(lines, s) {
			return "", nil
		}
		lines = append(lines, s)
		return s, nil
	}
	fns := make([]func(string) (string, error), T.NTerms())
	for i, _ := range fns {
		fns[i] = f
	}
	h := T.Headers()
	T.SetMany(h, fns)
}

// Returns the same function that will delete lines that are identical to the previous non-repeated
// line in the file. So, if 2 lines are repeated, but there are lines in between them, the repeated
// ones will _not_ be deleted. For that, see DelRepeatedFunctions. The reason for adding this function
// is that it requires less memory.
func DelDuplicatedFunctions(T TermSelect) {
	prevline := "2@@$)(/^*9HhL.Goは最高のプログラミング言語です" //just a line that is unlikely to be found in a topology file
	f := func(s string) (string, error) {
		if s == prevline {
			return "", nil
		}
		prevline = s
		return s, nil
	}
	fns := make([]func(string) (string, error), T.NTerms())
	for i, _ := range fns {
		fns[i] = f
	}
	h := T.Headers()
	T.SetMany(h, fns)
}

// The terms string must be already formatted in Gromacs format for each term.
// The slice of slices must contain one slice per header (which can be empty).
// in the alphabetical order of the headers.
// if you want to set the non-selected terms to 'do nothing' you can do so before
// calling this, with the method 'SetAll'
// and each slice contains a term to be added for that type.
func AddTermFunctions(T TermSelect, terms2add [][]string) {
	sel := T.Headers(true)
	slices.Sort(sel)
	for i, k := range sel {
		call := 0
		T.m[k] = func(string) (string, error) {
			if call == 0 {
				call++ //I _think_ each 'version' of the variable call will be 'closured' in each function
				//so each function will only work once.
				//I admit it's not very efficient, but I find it way clearer than dealing with it at
				//the function applier level.
				return strings.Join(terms2add[i], ""), nil
			}
			return "", nil
		}
	}
}

// Returns a TermSelect with funcions in atoms and vsites that delete all atoms and vsites
// with indexes present in todel. todel must be 1-based indexes. The functions do _not_
// adjust the indexes of the atoms not removed (see the ShiftAtomNumberFunctions)
// it also do _not_ remove bonded terms involving the atoms.
func DeleteAtomsFunction(todel []int) TermSelect {
	f := func(s string) (string, error) {
		at, err := strconv.Atoi(fi(s)[0])
		if err != nil {
			return "", fmt.Errorf("Can't parse numbers in line %s %w", s, err)
		}
		if slices.Contains(todel, at) {
			return "", nil
		}
		return s, nil
	}
	T := NewTermSelect()
	T.SetAll() //set everything to a 'do nothing' function
	//so everything we don't set explicily will just use that.
	T.Set("atoms", f)
	T.Set("vsitesn", f)

	return T
}

// Returns a TermSelect with functions that replace the term belonging to one of the headers in
// headers2edit and with atoms atoms (which must match in order, too, or in reverse order in
// 'reverse' is true), with the string newterm, unless the newterm string is "DONOTREPLACE" in which case
// the functions in the returned TermSelect do nothing. In addition EditTermFunction returns a function that, after
// the TermSelect is applied to a topology, will return a slice of strings with the matching terms.
func EditTermFunction(newterm string, atoms []int, headers2edit []string, reverse bool) (TermSelect, func() []string) {
	retterm := []string{}
	replace := func(s string) string {
		if newterm == "DONOTREPLACE" {
			return s
		}
		return newterm
	}
	adder := func() func(s string) (string, error) {
		return func(s string) (string, error) {
			st := fi(s)
			if len(atoms) <= 0 || len(atoms) > len(st) {
				return "", fmt.Errorf("0 or too many elements in atom list: %d", len(atoms))
			}
			///		fmt.Println(len(st), atoms, st, s) ////////////////////////////////
			ats, err := parseints(st[:len(atoms)]...)
			if err != nil {
				return "", fmt.Errorf("Can't parse numbers in line %s %w", s, err)
			}
			fmt.Println("ats, atoms", ats, atoms) ////////////////////////////////////////////
			if slices.Equal(ats, atoms) {
				retterm = append(retterm, s)
				return replace(s), nil
			}
			if reverse {
				slices.Reverse(ats)
				if slices.Equal(ats, atoms) {
					retterm = append(retterm, s)
					return replace(s), nil
				}

			}
			return s, nil
		}
	}
	T := NewTermSelect()
	h := T.Headers()
	T.SetAll()
	fmt.Println("header2edit", headers2edit) /////////////////////
	for _, v := range h {
		if slices.Contains(headers2edit, v) {
			T.Set(v, adder())
		}
	}
	return T, func() []string { return retterm }

}

// Returns a set of functions that will switch every atom index in tosub for the atom in the same position
// in replacement, in lines belonging to supported records of an itp file (as of now: atoms, bonds, constraints, angles, dihedrals,
// exclusions, virtual_sitesn.
// Some exceptional behavior:
// If len(tosub)==0 then replacement is expected to have 1 value, which will be added to all atoms.
// if len(tosub)==1 and that element is negative, then replacement is also expected to have one value, which will be
// added to all atoms with index larger than the absolute value of tosub[0]
// I admit this 'exceptionsl behavior' is pretty hack-y, but I didn't want to repeat code so much.
// if you need to add/subtract something to all atoms greater than something, it's easier to use ShiftAtomNumberFunctions
func SwitchAtomNumbersFunctions(tosub, replacement []int) TermSelect {
	//	rep := strings.Replace
	change := func(i int) int {
		if len(tosub) == 0 {
			return i + replacement[0]
		}
		index := slices.Index(tosub, i)
		if index == -1 {
			return i
		}
		return replacement[index]
	}
	after := 0
	if len(tosub) == 1 && tosub[0] < 0 {
		after = -1 * tosub[0]
		tosub = []int{}
	}

	if len(tosub) != 0 && len(tosub) != len(replacement) {
		panic(sf("tosub and replacement should have equal lenght unless len(tosub)==0, but len(tosub): %d len(replacement): %d", len(tosub), len(replacement)))
	}
	T := NewTermSelect()

	adder := func(atoms int) func(s string) (string, error) {
		return func(s string) (string, error) {
			st := fi(s)
			if atoms <= 0 || atoms > len(st) {
				atoms = len(st) //all fields are atoms
			}
			///		fmt.Println(len(st), atoms, st, s) ////////////////////////////////
			ats, err := parseints(st[:atoms]...)
			if err != nil {
				return "", fmt.Errorf("Can't parse numbers in line %s %w", s, err)
			}
			for i, v := range ats {
				if v < after {
					continue
				}
				at2 := change(v)
				st[i] = sf("%3d", at2)
			}
			return strings.Join(st, " "), nil
		}
	}

	vsitesn := func(s string) (string, error) {
		st := fi(s)
		ats, err := parseints(st...)
		if err != nil {
			return "", fmt.Errorf("Can't parse numbers in line %s %w", s, err)
		}
		if ats[0] > after {
			at2 := change(ats[0])
			st[0] = sf("%3d", at2)

		}
		for i, v := range ats[2:] {
			if v <= after {
				continue
			}
			at2 := change(v)
			st[i+2] = sf("%3d", at2)
		}
		return strings.Join(st, " "), nil
	}
	T.SetAll() //set everything to a 'do nothing' function
	//so everything we don't set explicily will just use that.
	T.Set("atoms", adder(1))
	T.Set("bonds", adder(2))
	T.Set("angles", adder(3))
	T.Set("constraints", T.Get("bonds"))
	T.Set("dihedrals", adder(4))
	T.Set("exclusions", adder(-1))
	T.Set("vsitesn", vsitesn)

	return T
}

// Returns a set of functions that will add or subtract toadd (depending on its sign) to every atom index greater
// than after, in lines belonging to supported records of an itp file (as of now: atoms, bonds, constraints, angles, dihedrals,
// exclusions, virtual_sitesn
func ShiftAtomNumbersFunctions(after, toadd int) TermSelect {
	return SwitchAtomNumbersFunctions([]int{-1 * after}, []int{toadd})
}

//The next one is not a utility function, as it is applied by itself, without the FuncApplier function, but I'll keep it in this file

// Merges the '[]' sections of two topologies. The sections in both must be in the same order.
// It adds all
func MergeTopologies(top1, top2 StringReader, target io.StringWriter) error {
	header := NewTopHeader()
	var fl string
	var err error
	next_header2 := ""
	for fl, err = top1.ReadString('\n'); err == nil; fl, err = top1.ReadString('\n') {
		ls := strings.Split(fl, ";") //remove comments
		l := ls[0]
		target.WriteString(fl)
		//	fmt.Println("casha", fl) /////////////
		//every time we reach a header in top1, we start looking for the header in top2
		//if we reach a header, or we have previously reached a header (next_header2),
		//and it's the same we read in top2, then we start inserting
		//whatever is in top2 after that header, until we encounter the next header.
		if header.Is(l) {
			writing2 := false
			h1 := header.Which(l)
			//	fmt.Println("headerq1", l, h1) ///////////////////
			if h1 == "moleculetype" {
				continue
			}
			var err2 error
			var fl2 string
			for fl2, err2 = top2.ReadString('\n'); err2 == nil; fl2, err2 = top2.ReadString('\n') {

				ls2 := strings.Split(fl2, ";") //remove comments
				l2 := ls2[0]
				if header.Is(l2) && header.Which(l2) == h1 {
					writing2 = true

					//			fmt.Println("headerq2", l2) ///////////////////
					continue
				}
				if strings.HasPrefix(l2, "#") {
					target.WriteString(fl2)
					continue
				}
				if next_header2 == h1 || writing2 {
					//			fmt.Println("next_header1ql", next_header2) ///////////////////

					if header.Is(l2) {
						next_header2 = header.Which(l2)
						writing2 = false
						break
					}

					//			fmt.Println("cashi", fl2) //////////
					target.WriteString(fl2)
				}
			}
			//		v, ok := top2.(*TopInMem)
			//		if ok {
			//			v.Reset()
			//		}
			if err2 != nil && !errors.Is(err2, io.EOF) {
				//		fmt.Println("this error shouldn't be here", err2) ///////////
				return err2
			}
		}
	}
	if !errors.Is(err, io.EOF) {
		return err
	}
	return nil
}
