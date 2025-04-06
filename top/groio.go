package top

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"os"
	"slices"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
)

type cond struct {
	reading bool
}

func newCond() *cond {
	c := new(cond)
	c.reading = true
	return c
}

// a function to read conditional parts of gromacs topologies
// depending on the defined flags that should be in 'defines'
func (c *cond) read(line string, defines []string) bool {
	if strings.HasPrefix(line, "#ifdef") {
		if slices.Contains(defines, fi(line)[0]) {
			c.reading = true
			return false
		} else {
			c.reading = false
			return false
		}
	}
	if strings.HasPrefix(line, "#else") {
		c.reading = !c.reading
		return false
	}
	if strings.HasPrefix(line, "#endif") {
		c.reading = true
		return false
	}
	return c.reading
}

//The high-level functions

// FillFF will fill the receiver with data from the given StringReader which must be
// in Gromacs itp/top format. If non bonded (Lennard-Jones) terms are present it
// will interpret them as sigma/epsilon if sigmaep is true, as C6/C12 otherwise
// if followIncludes values are given, the first one will determine whether #include
// statements will trigger opening and reading the included file(s).
// NOTE: Many Gromacs headers are supported, but not all. Notably, virtual_sitesX
// headers, where X is between 1-4 are currently _not_ supported.
// Only the generic virtual_sitesn are.
func (F *FF) Fill(r StringReader, followIncludes bool, defines ...string) error {
	follow := followIncludes
	var err error
	var s string
	read := newCond()
	h := NewTopHeader()
	if F.Mol == nil {
		return fmt.Errorf("Can't read topology if no molecule is loaded on the FF object")
	}
	for s, err = r.ReadString('\n'); err == nil; s, err = r.ReadString('\n') {
		s = cleanString(s)
		if s == "" {
			continue
		}
		if !read.read(s, defines) {
			continue
		}
		//This should allow us to 'follow' include files. Probably a risky thing to use.
		if strings.Contains(s, "#include") && follow {
			f := fi(s)
			fname := strings.Trim(f[len(f)-1], "\"'")
			file, err := os.Open(fname)
			if err != nil {
				return fmt.Errorf("Failed to include file: %s. Error: %w", file, err)
			}
			defer file.Close()
			reader := bufio.NewReader(file)
			err = F.Fill(reader, follow)
			if err != nil {
				return fmt.Errorf("Failed to include file: %s. Error: %w", file, err)
			}
			continue
		}
		if h.Is(s) {
			F.currentHeader = h.Which(s)
			continue
		}
		var T *Term
		var vs *VSite
		var att *AtomType
		var LJ *LJPair

		switch F.currentHeader {

		case "atoms":
			err = F.AtomDataFromGro(s)
		case "bonds":
			T, err = TermFromGro(s, F.currentHeader)
			F.Bonds = append(F.Bonds, T)
		case "constraints":
			T, err = TermFromGro(s, F.currentHeader)
			F.Constraints = append(F.Constraints, T)
		case "angles":
			T, err = TermFromGro(s, F.currentHeader)
			F.Angles = append(F.Angles, T)
		case "impropers":
			T, err = TermFromGro(s, F.currentHeader)
			F.Impropers = append(F.Impropers, T)
		case "dihedrals":
			T, err = TermFromGro(s, F.currentHeader)
			F.Dihedrals = append(F.Dihedrals, T)
		//note that I don't support other vsites right now.
		case "vsitesn":
			vs, err = VSitesNFromGro(s)
			F.VSites = append(F.VSites, vs)
		case "exclusions":
			err = F.ExclusionsFromGro(s)
		case "atomtypes":
			att, err = AtomTypeFromGro(s, F.SigmaEpsilon)
			F.ATypes = append(F.ATypes, att)
		case "nonbond":
			LJ, err = LJPairFromGro(s, F.SigmaEpsilon)
			F.LJ = append(F.LJ, LJ)
		case "defaults":
			f := fi(s)
			if f[0] == "1" && f[1] == "2" {
				F.SigmaEpsilon = true
			} else {
				F.SigmaEpsilon = false
			}
		default:
			continue

		}

		if err != nil {
			return fmt.Errorf("Couldn't read header %s. Line: %s. Error: %w", F.currentHeader, s, err)
		}
	}
	if errors.Is(err, io.EOF) {
		err = nil
	}
	return err
}

func (F *FF) AllToGro(r io.StringWriter) (err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()

	if len(F.ATypes) > 0 {
		_, err = r.WriteString("[ atomtypes ]\n")
		qerr(err)
		printGro(r, F.ATypes)
	}
	if len(F.LJ) > 0 {
		_, err = r.WriteString("\n[ nonbond_params ]\n")
		qerr(err)
		printGro(r, F.LJ)
	}
	_, err = r.WriteString("\n[ atoms ]\n")
	qerr(err)
	for i := 0; i < F.Mol.Len(); i++ {
		_, err = r.WriteString(F.Atom2Gro(i))
		qerr(err)
	}
	_, err = r.WriteString("\n[ bonds ]\n")
	printGro(r, F.Bonds)
	_, err = r.WriteString("\n[ constraints ]\n")
	printGro(r, F.Constraints)
	_, err = r.WriteString("\n[ exclusions ]\n")
	printExcluGro(r, F.Exclusions)
	_, err = r.WriteString("\n[ angles ]\n")
	printGro(r, F.Angles)
	_, err = r.WriteString("\n[ impropers ]\n")
	printGro(r, F.Impropers)
	_, err = r.WriteString("\n[ dihedrals ]\n")
	printGro(r, F.Dihedrals)
	_, err = r.WriteString("\n[ virtual_sitesn ]\n")
	//while the ToGro() signature for the interface
	//printGo takes requires returning an error, the only
	//implementation that actually has an error to return
	//is that of VSites.
	err = printGro(r, F.VSites)
	qerr(err)

	return nil
}

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

type groer interface {
	ToGro() (string, error)
}

func printGro[G ~[]E, E groer](r io.StringWriter, g G) error {
	for _, v := range g {
		m, e := v.ToGro()
		if e != nil {
			return e
		}
		_, e = r.WriteString(m)
		if e != nil {
			return e
		}
	}
	return nil
}

func printExcluGro(r io.StringWriter, g [][]int) error {
	for _, v := range g {
		v1 := exclusion(v)
		m, e := v1.ToGro()
		if e != nil {
			return e
		}
		_, e = r.WriteString(m)
		if e != nil {
			return e
		}
	}
	return nil
}

// Returns a string without gromacs comments (sequences starting with ';'),
// trailing and leading spaces, tabs and newlines
func cleanString(s string) string {
	f := strings.Split(s, ";")[0]
	return strings.Trim(f, "\n\t ")

}

func (F *FF) ExclusionsFromGro(s string) error {
	ex, err := parseints(fi(s)...)
	if err != nil {
		return err
	}
	F.Exclusions = append(F.Exclusions, ex)
	return nil

}

type exclusion []int

func (e exclusion) ToGro() (string, error) {
	ret := make([]string, 0, len(e)+1)
	for _, v := range e {
		ret = append(ret, sf("%4d", v))
	}
	ret = append(ret, "\n")
	return strings.Join(ret, " "), nil

}

// Adds the data in the gromacs-topology atom-section string to the atom with index index[0] in the
// molecule in ff. IF index is not given, the function will search the molecule to add the data to the
// atom that matches the ID on the topology string.
func (F *FF) AtomDataFromGro(s string, index ...int) (err error) {
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
func (F *FF) Atom2Gro(i int) string {
	A := F.Mol.Atom(i)
	return fmt.Sprintf("    %5d  %4s %5d  %4s  %4s %5d  %6.4f \n", A.ID, A.Symbol, A.MolID, A.MolName, A.Name, A.ID, A.Charge)
}

// Returns a term containing the information in the GromacsTop-formatted string s,
// given that the string is part of the header header.
func TermFromGro(s, header string) (T *Term, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	T = new(Term)
	l := strings.Fields(cleanString(s))
	T.OneBased = 0

	//first the 'special' cases.

	//r-b terms. don't have equilibrium values or force constants.
	if header == "dihedrals" && len(l) == 11 {
		// Ryckaert-Bellemans entonces
		T.IDs, err = parseints(l[:4]...)
		qerr(err)
		T.FuncType = 3
		T.RB, err = parsefloats(l[5:]...)
		qerr(err)
		if len(T.RB) != 6 {
			err = fmt.Errorf("R-B term detected but read %d parameters instead of the 6 expected", len(T.RB))
			T = nil
			return
		}
		return
	}
	//constraints are especial because they lack force constants.
	if header == "constraints" {
		T.IDs, err = parseints(l[:2]...)
		qerr(err)
		T.Constraint = true
		var ft int
		ft, err = strconv.Atoi(l[2])
		qerr(err)
		T.FuncType = uint(ft)
		T.Eq, err = strconv.ParseFloat(l[3], 64) //could be 3 if there is a function type but I don't think there is. Check
		qerr(err)
		return
	}
	//everything else is: atom1 ... atom(ats) functype eqval fconst
	var ats int = -1
	switch header {
	case "bonds":
		ats = 2
	case "angles":
		ats = 3
	case "impropers":
		ats = 4
	case "dihedrals":
		ats = 4
	}

	T.IDs, err = parseints(l[:ats]...)
	T.Constraint = false
	var ft int
	ft, err = strconv.Atoi(l[ats])
	qerr(err)
	T.FuncType = uint(ft)
	T.Eq, err = strconv.ParseFloat(l[ats+1], 64) //could be 3 if there is a function type but I don't think there is. Check
	qerr(err)
	T.K, err = strconv.ParseFloat(l[ats+2], 64) //could be 3 if there is a function type but I don't think there is. Check
	qerr(err)
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
func (T *Term) ToGro() (string, error) {
	ret := make([]string, 0, 15)
	ret = append(ret, T.writeAtoms())
	if T.Vsite {
		return "", nil //placeholder
	}
	if !T.Constraint {
		ret = append(ret, fmt.Sprintf("%1d", T.FuncType))
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
	ret = append(ret, "\n")

	return strings.Join(ret, " "), nil
}

// Returns a *VSite witht he information in
// the string s containing a virtual_sitesn line
func VSitesNFromGro(s string) (vsit *VSite, err error) {
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
		for _, v := range f[2:] {
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
	ret := make([]string, 0, 8)
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

// Reads a string with the appropriate gromacs topology format
// to return a pointer to AtomType. If sigmaep is true, it transforms
// the data in the string from sigma/epsilon to c6/c12
func AtomTypeFromGro(s string, sigmaep bool) (ret *AtomType, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("Couldn't read atom type from string. Error: %s String:%s", r, s)
		}
	}()
	s = cleanString(s)
	f := fi(s)
	ret = new(AtomType)
	ret.Name = f[0]
	ret.Mass, err = strconv.ParseFloat(f[1], 64)
	qerr(err)
	ret.Charge, err = strconv.ParseFloat(f[2], 64)
	qerr(err)
	ret.Ptype = f[3]

	ret.C6, ret.C12, err = c6c12OrSigmaEpsilon(f[4], f[5], sigmaep)
	ret.SigmaEpsilon = sigmaep
	qerr(err)
	return ret, err
}

func (A *AtomType) ToGro() (string, error) {
	var c6, c12 float64
	if A.SigmaEpsilon {
		c6, c12 = c6c12ToSigmaepsilon(A.C6, A.C12)

	} else {
		c6, c12 = A.C6, A.C12
	}
	return sf("%5s %5.3f %5.3f %1s %6.4e %6.4ef", A.Name, A.Mass, A.Charge, A.Ptype, c6, c12), nil
}

func LJPairFromGro(s string, sigmaep bool) (ret *LJPair, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("Couldn't read atom type from string. Error: %s String:%s", r, s)
		}
	}()
	s = cleanString(s)
	f := fi(s)
	ret = new(LJPair)
	ret.Names[0] = f[0]
	ret.Names[1] = f[1]
	ret.FuncType, err = strconv.Atoi(f[2])
	qerr(err)
	ret.C6, ret.C12, err = c6c12OrSigmaEpsilon(f[4], f[5], sigmaep)
	qerr(err)
	ret.SigmaEpsilon = sigmaep
	return ret, err

}

func (L *LJPair) ToGro() (string, error) {
	var c6, c12 float64
	if L.SigmaEpsilon {
		c6, c12 = c6c12ToSigmaepsilon(L.C6, L.C12)

	} else {
		c6, c12 = L.C6, L.C12
	}
	return sf("%5s %5s %1d %6.4e %6.4e\n", L.Names[0], L.Names[1], L.FuncType, c6, c12), nil
}

func c6c12OrSigmaEpsilon(num1, num2 string, sigmaepsilon bool) (float64, float64, error) {
	var c6, c12 float64
	var err error
	c6, err = strconv.ParseFloat(num1, 64)
	if err != nil {
		return -1, -1, err
	}
	c12, err = strconv.ParseFloat(num2, 64)
	if err != nil {
		return -1, -1, err
	}
	if sigmaepsilon {
		c6, c12 = sigmaepsilonToc6c2(c6, c12)

	}
	return c6, c12, nil

}
