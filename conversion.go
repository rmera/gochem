/*
 * conversion.go, part of gochem.
 *
 *
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package chem

//This provides useful conversion factors and other constants

//Conversions
const (
	Deg2Rad = 0.0174533
	Rad2Deg = 1 / 0.0174533
	H2Kcal  = 627.509 //HArtree 2 Kcal/mol
	Kcal2H  = 1 / 627.509
	KJ2Kcal = 1 / 4.184
	Kcal2KJ = 4.184
	A2Bohr  = 1.889725989
	Bohr2A  = 1 / 1.889725989
	EV2Kcal = 23.061
	Kcal2EV = 1 / 23.061
)

//Others
const (
	CHDist = 1.098               //C(sp3)--H distance in A
	KB     = 1.380649e-26        // Boltzmann constant kJ/K
	R      = 8.31446261815324e-3 // kJ/(K*mol)
	NA     = 6.02214076e+23
)
