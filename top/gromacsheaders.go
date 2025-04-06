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

package top

import (
	"fmt"
	"os/exec"
	"regexp"
	"strconv"
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
		"atomtypes":    regexp.MustCompile(`\[\p{Zs}atomtypes\p{Zs}*\]`),
		"nonbond":      regexp.MustCompile(`\[\p{Zs}nonbond_params\p{Zs}*\]`),
		"defaults":     regexp.MustCompile(`\[\p{Zs}*defaults\p{Zs}*\]`),
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
