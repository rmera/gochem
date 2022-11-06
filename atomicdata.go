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
 * goChem is currently developed at the Universidad de Santiago de Chile
 * (USACH)
 *
 */

package chem

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	"strconv"
	"strings"
)

//A map for assigning mass to elements.
//Note that just common "bio-elements" are present
var symbolMass = map[string]float64{
	"H":  1.0,
	"C":  12.01,
	"O":  16.00,
	"N":  14.01,
	"P":  30.97,
	"S":  32.06,
	"Se": 78.96,
	"K":  39.1,
	"Ca": 40.08,
	"Mg": 24.30,
	"Cl": 35.45,
	"Na": 22.99,
	"Cu": 63.55,
	"Zn": 65.38,
	"Co": 58.93,
	"Fe": 55.84,
	"Mn": 54.94,
	"Cr": 51.996,
	"Si": 28.08,
	"Be": 9.012,
	"F":  18.998,
	"Br": 79.904,
	"I":  126.90,
}

//A map for assigning covalent radii to elements
//Values from Cordero et al., 2008 (DOI:10.1039/B801115J)
//Note that just common "bio-elements" are present
var symbolCovrad = map[string]float64{
	"H":  0.4,  // 0.31 I altered this one. Since H always has only one bond, it doesn't matter if I set a longer radius, the extra bonds will get eliminated later.
	"C":  0.76, //the sp3 radius
	"O":  0.66,
	"N":  0.71,
	"P":  1.07,
	"S":  1.05,
	"Se": 1.2,
	"K":  2.03,
	"Ca": 1.76,
	"Mg": 1.41,
	"Cl": 1.02,
	"Na": 1.66,
	"Cu": 1.32,
	"Zn": 1.22,
	"Co": 1.5,  // hs
	"Fe": 1.52, //hs
	"Mn": 1.61, //hs
	"Cr": 1.39,
	"Si": 1.11,
	"Be": 0.96,
	"F":  0.57,
	"Br": 1.2,
	"I":  1.39,
}

//A map for assigning van der Waals radii to elements
//Values from 10.1021/j100785a001 and 10.1021/jp8111556
//metal radii from 10.1023/A:1011625728803
//Note that just common "bio-elements" are present
var symbolVdwrad = map[string]float64{
	"H":  1.10, // 0.31 I altered this one. Since H always has only one bond, it doesn't matter if I set a longer radius, the extra bonds will get eliminated later.
	"C":  1.70, //the sp3 radius
	"O":  1.52,
	"N":  1.55,
	"P":  1.80,
	"S":  1.80,
	"Se": 1.90,
	"K":  2.75,
	"Ca": 2.31,
	"Mg": 1.73,
	"Cl": 1.75,
	"Na": 2.27,
	"Cu": 2.00,
	"Zn": 2.02,
	"Co": 1.95,
	"Fe": 1.96,
	"Mn": 1.96,
	"Cr": 1.97,
	"Si": 2.10,
	"Be": 1.53,
	"F":  1.47,
	"Br": 1.83,
	"I":  1.98,
}

//A map for checking that atoms don't
//have too many bonds. A value of 0 means
//undefined, i.e. that this atom shouldn't
//be checked for max bonds. I decided not to define it
var symbolMaxBonds = map[string]int{
	"H":  1, //this is the only one truly important.
	"C":  4,
	"O":  2,
	"N":  0, //undefined
	"P":  0,
	"S":  0,
	"Se": 0,
	"Be": 0,
	"F":  1,
	"Br": 1,
	"I":  1,
}

//ALL the Following is _not_ tested.

type RA2VdW struct {
	a   *regexp.Regexp
	vdw float64
}

//Martini 3 bible, p 42
var Martini3TypesvdW []*RA2VdW = []*RA2VdW{&RA2VdW{a: regexp.MustCompile("^S"), vdw: 2.30}, &RA2VdW{a: regexp.MustCompile("^T"), vdw: 1.91}, &RA2VdW{a: regexp.MustCompile("^(Q|P|N|C|X|D)"), vdw: 2.64}, &RA2VdW{a: regexp.MustCompile("W"), vdw: 1.4 * 4}}

//Martini 3 bible, p 42
var Martini3VdW []*RA2VdW = []*RA2VdW{&RA2VdW{a: regexp.MustCompile(".*BB"), vdw: 2.30}, &RA2VdW{a: regexp.MustCompile("(CYS|ALA|SER|HIS|HIH|TYR|TRP)SC"), vdw: 1.91}, &RA2VdW{a: regexp.MustCompile("PHESC1"), vdw: 2.3}, &RA2VdW{a: regexp.MustCompile("PHESC(2|3)"), vdw: 1.91}, &RA2VdW{a: regexp.MustCompile("W"), vdw: 1.4 * 4}}

//Where REGEX is a regular expression to be matched with atom residue and atomname (example: TYRCA)  and radius is the
//corresponding vdW radius for the matching atoms, in A. The first REGEX matched will be used
//for each atom, so the order of the lines in the file might be significant.
func VdwResAndNAmeFromFile(fname string, mol Atomer) ([]*RA2VdW, error) {
	f, err := os.Open(fname)
	if err != nil {
		return nil, err
	}
	buf := bufio.NewReader(f)

	radii := make([]*RA2VdW, 0, mol.Len())
	var line string
	for line, err = buf.ReadString('\n'); err == nil; line, err = buf.ReadString('\n') {
		f := strings.Fields(line)
		if len(f) < 2 {
			err = &CError{msg: fmt.Sprintf("vdW file %s improperly formatted", fname), deco: []string{}}
			break
		}
		re, err2 := regexp.Compile(f[0])
		if err2 != nil {
			err = err2
			break
		}
		val, err3 := strconv.ParseFloat(f[1], 64)
		if err3 != nil {
			err = err3
			break
		}
		radii = append(radii, &RA2VdW{a: re, vdw: val})
	}
	if err != nil && err.Error() != "EOF" {
		return nil, err
	}

	for i := 0; i < mol.Len(); i++ {
		v := mol.Atom(i)
		for _, w := range radii {
			if w.a.MatchString(v.Name) {
				v.Vdw = w.vdw
				break
			}
		}
	}
	return radii, nil
}
