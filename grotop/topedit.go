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
	"bufio"
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
		"atoms":        regexp.MustCompile(`\[\p{Zs}*atoms\p{Zs}*\]`),
		"constraints":  regexp.MustCompile(`\[\p{Zs}*constraints\p{Zs}*\]`),
		"bonds":        regexp.MustCompile(`\[\p{Zs}*bonds\p{Zs}*\]`),
		"angles":       regexp.MustCompile(`\[\p{Zs}*angles\p{Zs}*\]`),
		"vsites1":      regexp.MustCompile(`\[\p{Zs}*virtual_sites1\p{Zs}*\]`),
		"vsites2":      regexp.MustCompile(`\[\p{Zs}*virtual_sites2\p{Zs}*\]`),
		"vsites3":      regexp.MustCompile(`\[\p{Zs}*virtual_sites3\p{Zs}*\]`),
		"vsitesn":      regexp.MustCompile(`\[\p{Zs}*virtual_sitesn\p{Zs}*\]`),
		"exclusions":   regexp.MustCompile(`\[\p{Zs}*exclusions\p{Zs}*\]`),
		"molecules":    regexp.MustCompile(`\[\p{Zs}*molecules\p{Zs}*\]`),
		"dihedrals":    regexp.MustCompile(`\[\p{Zs}*dihedrals\p{Zs}*\]`),
		"moleculetype": regexp.MustCompile(`\[\p{Zs}moleculetype\p{Zs}*\]`),
	}
	T.set = true

}

func (T *topHeader) delcomments(line string) string {
	return cleanString(line)
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

func (T TermSelect) GetErr(header string) (func(string) (string, error), error) {
	val, ok := T.m[header]
	if !ok {
		return nil, fmt.Errorf("TermSelect.Get: Attempted to get unsuported value: %s", header)
	}
	if val == nil {
		return val, fmt.Errorf("Warning: Function set to nil")
	}

	return val, nil
}

// gets the function set to a header. Panics if the heeader is not present
func (T TermSelect) Get(header string) func(string) (string, error) {
	val, ok := T.m[header]
	if !ok {
		panic(fmt.Sprintf("grotop/TermSelect/MustGet: Attempted to get unsuported header: %s", header))
	}
	return val
}

func (T TermSelect) GetOrDefault(header string) func(string) (string, error) {
	val, ok := T.m[header]
	if !ok {
		panic(fmt.Sprintf("grotop/TermSelect/MustGet: Attempted to get unsuported header: %s", header))
	}
	if val == nil {
		return func(s string) (string, error) { return s, nil }
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

// Set all headers to a do-nothing function. Only the first given functionw
// will be used to set all headers. If nothing is given, a 'do nothing' function
// which returns the same string given, will be used.
func (t TermSelect) SetAll(f ...func(string) (string, error)) {
	fu := func(s string) (string, error) { return s, nil }
	if len(f) > 0 && f[0] != nil {
		fu = f[0]
	}
	for k, _ := range t.m {
		t.m[k] = fu
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

func TopInMemFromFile(fname string) (*TopInMem, error) {
	T := new(TopInMem)
	T.t = make([]string, 0, 10)
	f, err := os.Open(fname)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	re := bufio.NewReader(f)
	var l string
	for l, err = re.ReadString('\n'); err == nil; l, err = re.ReadString('\n') {
		ll := strings.TrimSuffix(l, "\n")
		T.t = append(T.t, ll)
	}
	if err != nil && errors.Is(err, io.EOF) {
		err = nil
	}
	//	fmt.Println(T.t) ////////////////////////////////////
	return T, err
}

// Returns a deep copy of the topology
func (t *TopInMem) Copy() *TopInMem {
	s := make([]string, len(t.t))
	copy(s, t.t)
	return NewTopInMem(s)
}

// Resets the reader to start from the first line
func (t *TopInMem) Reset() {
	t.i = 0
}

// Resets the reader to start from the first line
func (t *TopInMem) Len() int {
	return len(t.t)
}

// Adds a string to the topology
func (t *TopInMem) WriteString(s string) (int, error) {
	t.t = append(t.t, s)
	return len(s), nil

}

// Replaces the last-read string in the topology for s
func (t *TopInMem) ReplaceString(s string) {
	t.t[t.i-1] = s
}

// Returns the next line in the topology. Note that the byte argument is
// not used, you can't choose how much you want to read, it's always the
// full next line (unlike in the bufio.Reader ReadString method).
func (t *TopInMem) ReadString(byte) (string, error) {
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
		_, err = f.WriteString(v + "\n")
		if err != nil {
			return fmt.Errorf("grotop/TopInMem.WriteToFile: Couldn't write %d-th line to file: %w", i+1, err)
		}
	}
	return nil
}

// Removes unneeded whitelines (including lines containing only a ';') from t.
// it allocates a new []string of the same size as the original to do the change
// so not very memory efficient.
func (t *TopInMem) Clean() {
	newt := make([]string, 0, len(t.t))
	header := NewTopHeader() //will probably move thise to the structure
	//so it's not re-created every time you call this function.
	for i, v := range t.t {
		if strings.Trim(v, "\t; \n") == "" { //if the string is nothing but white spaces and linejumps
			if i < len(t.t)-1 && !header.Is(t.t[i+1]) {
				continue
			}
		}
		newt = append(newt, v)
	}
	t.t = newt
}

type StringReader interface {
	ReadString(byte) (string, error)
}
type StringReplacer interface {
	ReplaceString(string)
}
type StringReaderReplacer interface {
	StringReader
	StringReplacer
}

// Transforms each supported section of an Gromacs itp/top file with corresponding function in the
// given TermSelect. Writes the modified trajectory to a StringWriter, and returns an error or nil.
// if a function in the TermSelect map is nil, FuncApplier will apply a 'do nothing' function, that
// returns the same string given, instead.
func FuncApplier(top StringReaderReplacer, T TermSelect) error {
	header := NewTopHeader()
	currentfunc := T.GetOrDefault("default")
	currentheader := ""
	var fl string
	var err error
	for fl, err = top.ReadString('\n'); err == nil; fl, err = top.ReadString('\n') {
		if strings.HasPrefix(fl, ";") {
			top.ReplaceString(fl)
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
			li, err := T.GetOrDefault("nonskiplines")(l)
			if err != nil {
				return err
			}
			top.ReplaceString(li + co)
			continue
		}
		if strings.Contains(l, "#") {
			top.ReplaceString(l + co)
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

			top.ReplaceString(l + co)
			continue
		}
		li, err := currentfunc(l)
		if err != nil {
			return fmt.Errorf("Problem with header %s: %w", currentheader, err)
		}
		top.ReplaceString(li + co)

	}
	if errors.Is(err, io.EOF) {
		err = nil
	}
	return nil
}
