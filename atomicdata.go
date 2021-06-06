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
	"H":  0.31,
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

//A map for checking that atoms don't
//have too many bonds. A value of 0 means
//undefined, i.e. that this atom shouldn't
//be checked for max bonds.
var symbolMaxBonds = map[string]int{
	"H":  1,
	"C":  4, //the sp3 radius
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