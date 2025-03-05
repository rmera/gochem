/*
 * main.go, part of goChem
 *
 *
 * Copyright 2024 Raul Mera  <rmeraa{at}academicos(dot)uta(dot)cl>
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *
 */

/*
Set of functions to automate some simple editing of top files
*/

package gro

import (
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"regexp"
	"slices"
	"strconv"
	"strings"
)

// Utility functions

func qerr(err error) {
	if err != nil {
		panic(err.Error())
	}
}

func delcomment(line string) string {
	if strings.HasPrefix(line, ";") {
		return ""
	}
	return strings.Split(line, ";")[0] //remove comments
}

// yeah, yeah, it's ugly. If it makes you feel better, it's supposed to be temporary,
// while I'm still figuring out which commands do I need.
func runcq(command string, a ...interface{}) {
	w := exec.Command("sh", "-c", fmt.Sprintf(command, a...))
	w.Run()
}

func parseints(s ...string) ([]int, error) {
	r := make([]int, 0, len(s))
	for _, v := range s {
		i, err := strconv.Atoi(v)
		if err != nil {
			return nil, err
		}
		r = append(r, i)
	}
	return r, nil
}

func parsefloats(s ...string) ([]float64, error) {
	r := make([]float64, 0, len(s))
	for _, v := range s {
		i, err := strconv.ParseFloat(v, 64)
		if err != nil {
			return nil, err
		}
		r = append(r, i)
	}
	return r, nil
}

type Atom struct {
	ID      int
	MolID   int
	Name    string
	PDBName string
	MolName string
	Charge  float64
}

func ReadAtom(s string) (err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	A := new(Atom)
	l := strings.Fields(delcomment(s))
	A.ID, err = strconv.Atoi(l[0])
	qerr(err)
	A.Name = l[1]
	A.MolID, err = strconv.Atoi(l[2])
	qerr(err)
	A.MolName = l[3]
	A.PDBName = l[4]
	A.Charge, err = strconv.ParseFloat(l[6], 64)
	qerr(err)
	return nil
}

func (A *Atom) String() string {
	return fmt.Sprintf("    %5d  %4s %5d  %4s  %4s %5d  %6.4f \n", A.ID, A.Name, A.MolID, A.MolName, A.PDBName, A.ID, A.Charge)
}

type Term struct {
	Functype   uint
	Atoms      []int
	OneBased   bool
	K          float64
	Eq         float64
	Constraint bool
	Vsite      bool
	RB         []float64
}

func ReadTerm(s string) (T *Term, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	T = new(Term)
	l := strings.Fields(delcomment(s))
	if len(l) == 11 {
		// Ryckaert-Bellemans entonces
		T.Atoms, err = parseints(l[:4]...)
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
	if len(l) == 3 {
		//constraints
		T.Atoms, err = parseints(l[:2]...)
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
	add := 0
	if !T.OneBased {
		add++
	}
	r := make([]string, 0, len(T.Atoms))
	for _, v := range T.Atoms {
		r = append(r, fmt.Sprintf("%4d", v+add))
	}
	return strings.Join(r, " ")

}

func (T *Term) String() string {
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

//The following is to be removed and placed in a separate "utils" module

type topHeader struct {
	wany      *regexp.Regexp
	vsitesany *regexp.Regexp
	spec      map[string]*regexp.Regexp
	set       bool
	line      string
}

func NewTopHeader() *topHeader {
	R := new(topHeader)
	R.Set()
	return R

}

func (T *topHeader) Set() {
	T.wany = regexp.MustCompile(`\[\p{Zs}*.*\p{Zs}*\]`)
	T.vsitesany = regexp.MustCompile(`\[\p{Zs}*virtual_sites[01234n]?\p{Zs}*\]`)
	T.spec = map[string]*regexp.Regexp{
		"atoms":       regexp.MustCompile(`\[\p{Zs}*atoms\p{Zs}*\]`),
		"constraints": regexp.MustCompile(`\[\p{Zs}*constraints\p{Zs}*\]`),
		"bonds":       regexp.MustCompile(`\[\p{Zs}*bonds\p{Zs}*\]`),
		"angles":      regexp.MustCompile(`\[\p{Zs}*angles\p{Zs}*\]`),
		"vsites1":     regexp.MustCompile(`\[\p{Zs}*virtual_sites1\p{Zs}*\]`),
		"vsites2":     regexp.MustCompile(`\[\p{Zs}*virtual_sites2\p{Zs}*\]`),
		"vsites3":     regexp.MustCompile(`\[\p{Zs}*virtual_sites3\p{Zs}*\]`),
		"vsitesn":     regexp.MustCompile(`\[\p{Zs}*virtual_sitesn\p{Zs}*\]`),
		"exclusions":  regexp.MustCompile(`\[\p{Zs}*exclusions\p{Zs}*\]`),
		"molecules":   regexp.MustCompile(`\[\p{Zs}*molecules\p{Zs}*\]`),
		"dihedrals":   regexp.MustCompile(`\[\p{Zs}*dihedrals\p{Zs}*\]`),
	}
	T.set = true

}

func (T *topHeader) delcomments(line string) string {
	return delcomment(line)
}

// Returns true if the line is a Gromacs header. It discards comments.
func (T *topHeader) Is(line string) bool {
	line = T.delcomments(line)
	return T.wany.MatchString(line)
}

func (T *topHeader) IsVirtualSites(line string) bool {
	line = T.delcomments(line)
	return T.vsitesany.MatchString(line)
}

func (T *topHeader) IsDihedrals(line string) bool {
	line = T.delcomments(line)
	return T.spec["dihedrals"].MatchString(line)
}

// Returns a string indicating which Gromacs top file header
// the line is, or an empty string if the line is not a header.
func (T *topHeader) Which(line string) string {
	line = T.delcomments(line)
	if !T.wany.MatchString(line) {
		return ""
	}
	for k, v := range T.spec {
		if v.MatchString(line) {
			return k
		}
	}
	return ""
}

type TermSelect struct {
	m map[string]func(string) (string, error)
}

func NewTermSelect() TermSelect {
	m := make(map[string]func(string) (string, error))
	m["atoms"] = nil
	m["bonds"] = nil
	m["angles"] = nil
	m["constraints"] = nil
	m["dihedrals"] = nil
	m["exclusions"] = nil
	m["vsitesn"] = nil
	m["nonskiplines"] = nil
	m["default"] = nil
	return TermSelect{m: m}
}

func (T TermSelect) NTerms() int {
	return len(T.m)
}

func (T TermSelect) GetErr(s string) (func(string) (string, error), error) {
	val, ok := T.m[s]
	if !ok {
		return nil, fmt.Errorf("TermSelect.Get: Attempted to get unsuported value: %s", s)
	}

	return val, nil
}
func (T TermSelect) Get(s string) func(string) (string, error) {
	val, ok := T.m[s]
	if !ok {
		panic(fmt.Sprintf("grotop/TermSelect/MustGet: Attempted to get unsuported value: %s", s))
	}
	return val
}

// Set each of the given headers to the corresponding string slice in s (headers and s must
// be of the same lenght. Note that you can set
// headers to 'nil', thus marking them as 'unselected'. Panics if header doesn't exist.
func (T TermSelect) Set(header string, s func(string) (string, error)) {
	_, ok := T.m[header]
	if !ok {
		panic(fmt.Sprintf("term %s not supported", s))
	}
	T.m[header] = s
}

// Same as Set but returns an error if header doesn't exist.
func (T TermSelect) SetErr(header string, s func(string) (string, error)) (err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("%s", r)
		}
	}()
	T.Set(header, s)
	return
}

// Set each of the given headers to the corresponding string slice in s (headers and s must
// be of the same lenght. Note that you can set
// headers to 'nil', thus marking them as 'unselected'.
func (T TermSelect) SetMany(headers []string, s []func(string) (string, error)) error {
	if headers == nil {
		return fmt.Errorf("TermSelect.Set: No term given")
	}
	if len(headers) > len(s) {
		return fmt.Errorf("TermSelect.Set: Not enough strings for the headers to be set")
	}
	var err error
	for i, v := range headers {
		_, ok := T.m[v]
		if !ok {
			if err == nil {
				err = fmt.Errorf("term %s not supported", v)
			} else {
				err = fmt.Errorf("term %s not supported - %s", v, err.Error())
			}
		}
		T.m[v] = s[i]

	}
	return err
}

// Set all headers to a do-nothing function
func (t TermSelect) SelectAll() {
	for k, _ := range t.m {
		t.m[k] = func(s string) (string, error) { return s, nil }
	}

}

// Returns all the headers sorted. If onlySelected is given and true,
// only the headers that have associated a non-nil function are returned.
func (T TermSelect) Headers(onlySelected ...bool) []string {
	var sel bool
	if len(onlySelected) > 0 {
		sel = onlySelected[0]
	}
	ret := make([]string, 0, 2)
	for k, v := range T.m {
		if !sel || (v != nil) {
			ret = append(ret, k)
		}
	}
	slices.Sort(ret)
	return ret
}

// Represents a topology stored in memory, as opposed to in a file.
type TopInMem struct {
	t []string
	i int
}

// Returns a new TopInMem, with the topology
// represented by the given slice of strings (each
// string must correspond to one line of the file, including
// the respective '\n'.
func NewTopInMem(t []string) *TopInMem {
	return &TopInMem{t: t, i: 0}
}

// Adds a string to the topology
func (t *TopInMem) WriteString(s string) (int, error) {
	t.t = append(t.t, s)
	return len(s), nil

}

// Returns the next line in the topology. Note that the byte argument is
// not used, you can't choose how much you want to read, it's always the
// full next line (unlike in the bufio.Reader ReadString method).
func (t *TopInMem) Readstring(byte) (string, error) {
	if t.i >= len(t.t) {
		t.i = 0 //you can re-start reading it.
		return "", io.EOF
	}
	t.i++
	return t.t[t.i-1], nil
}

func (t *TopInMem) WriteToFile(name string) error {
	f, err := os.Create(name)
	if err != nil {
		return fmt.Errorf("grotop/TopInMem.WriteToFile: %w", err)
	}
	for i, v := range t.t {
		_, err = f.WriteString(v)
		if err != nil {
			return fmt.Errorf("grotop/TopInMem.WriteToFile: Couldn't write %d-th line to file: %w", i+1, err)
		}
	}
	return nil
}

/************************************************
*
*      The actual utility functions
*
*************************************************/

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

// Returns a set of functions that will add (if add is true) or subtract toadd to every atom index greater
// than after, in lines belonging to supported records of an itp file (as of now: atoms, bonds, constraints, angles, dihedrals,
// exclusions, virtual_sitesn
func AddOrDelAtomFunctions(add bool, after, toadd int) TermSelect {
	fi := strings.Fields
	sf := fmt.Sprintf
	//	rep := strings.Replace
	change := func(i int) int {
		if add {
			return i + toadd
		} else {
			return i - toadd
		}
	}
	T := NewTermSelect()
	donothing := func(s string) (string, error) { return s, nil }

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
				if v <= after {
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

	T.Set("atoms", adder(1))
	T.Set("bonds", adder(2))
	T.Set("angles", adder(3))
	T.Set("constraints", T.Get("bonds"))
	T.Set("dihedrals", adder(4))
	T.Set("exclusions", adder(-1))
	T.Set("vsitesn", vsitesn)
	T.Set("nonskiplines", donothing)
	T.Set("default", donothing)

	return T
}

// Merges the '[]' sections of two topologies. The sections in both must be in the same order.
// It adds all
func MergeTopologies(top1, top2 StringReader, target io.StringWriter) error {
	header := NewTopHeader()
	var fl string
	var err error
	for fl, err = top1.ReadString('\n'); err == nil; fl, err = top1.ReadString('\n') {
		writing2 := false
		ls := strings.Split(fl, ";") //remove comments
		l := ls[0]
		target.WriteString(fl)
		//every time we reach a header in top1, we start looking for the header in top2
		//if we reach a header, or we have previously reached a header (next_header),
		//and it's the same we read in top2, then we start inserting
		//whatever is in top2 after that header, until we encounter the next header.
		next_header2 := ""
		if header.Is(l) {
			h1 := header.Which(l)
			var err2 error
			var fl2 string
			for fl2, err2 = top2.ReadString('\n'); err2 == nil; fl, err2 = top2.ReadString('\n') {
				ls2 := strings.Split(fl2, ";") //remove comments
				l2 := ls2[0]
				if header.Is(l2) && header.Which(l2) == h1 {
					writing2 = true
					continue
				}
				if strings.HasPrefix(l2, "#") {
					target.WriteString(fl2)
					continue
				}
				if next_header2 == h1 || writing2 {
					if header.Is(l2) {
						next_header2 = header.Which(l2)
						writing2 = false
						break
					}
					target.WriteString(fl2)
				}
			}
			if !errors.Is(err2, io.EOF) {
				return err2
			}
		}
	}
	if !errors.Is(err, io.EOF) {
		return err
	}
	return nil
}

type StringReader interface {
	ReadString(byte) (string, error)
}

// Transforms each supported section of an Gromacs itp/top file with corresponding function in the
// given map. Writes the modified trajectory to a StringWriter, and returns an error or nil.
func FuncApplier(top StringReader, target io.StringWriter, T TermSelect) error {
	header := NewTopHeader()
	currentfunc := T.Get("default")
	currentheader := ""
	var fl string
	var err error
	for fl, err = top.ReadString('\n'); err == nil; fl, err = top.ReadString('\n') {
		if strings.HasPrefix(fl, ";") {
			_, err = target.WriteString(fl)
			if err != nil {
				return err
			}
			continue
		}
		ls := strings.Split(fl, ";") //remove comments
		l := ls[0]
		co := ""
		if len(ls) > 1 {
			co = "; " + strings.Join(ls[1:], " ")
			//	co = strings.Replace(co, "\n", "", 1)
		}
		if len(strings.Fields(l)) == 0 {
			li, err := T.Get("nonskiplines")(l)
			if err != nil {
				return err
			}
			_, err = target.WriteString(li + co)
			if err != nil {
				return err
			}

			continue
		}
		if strings.Contains(l, "#") {
			_, err = target.WriteString(l + co)
			if err != nil {
				return err
			}

			continue
		}
		if header.Is(l) {
			fun, e := T.GetErr(header.Which(l))
			if e == nil {
				currentfunc = fun
				currentheader = header.Which(l)
				//return "", fmt.Errorf("couldn't find function %s in map", l)
			}

			_, err = target.WriteString(l + co)
			if err != nil {
				return err
			}

			continue
		}
		li, err := currentfunc(l)
		if err != nil {
			return fmt.Errorf("Problem with header %s: %w", currentheader, err)
		}
		_, err = target.WriteString(li + co)
		if err != nil {
			return err
		}

	}
	if errors.Is(err, io.EOF) {
		err = nil
	}
	return nil

}
