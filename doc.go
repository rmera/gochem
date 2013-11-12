/*
 * doc.go, part of gochem.
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

/*Package chem provides atom and molecule structures, facilities for reading and writing some
files used in computational chemistry and some functions for geometric manipulations and shape
indicators. 

Subdirectories of chem provide trajectory file reading, interaction with quantum chemistry programs
and Ramachandran plot capabilities.
chem Implements its own matrix type for coordinates, VecMatrix,based in github.com/gonum/matrix. 

Currently, each row of a VecMatrix represents one point in space. As this could change if 
github.com/gonum/matrix changes, we recomend prefering the Vec* methods over the Row* methods
when manipulating a VecMatrix. Vec* methods will change from row to column following the upstream
Gonum library.*/
package chem
