/*
 * doc.go, part of gochem.
 *
 * Copyright 2015 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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

/*Package v3 implements a Matrix type representing a row-major 3D matrix (i.e. a Nx3 matrix).
The v3.Matrix is used to represent the cartesian coordinates of sets of atoms in goChem.
It is based int gonum's (github.com/gonum) Dense type, with some additional restrictions
because of the fixed number of columns and with some additional functions that were found
useful for the purposes of goChem.

*/
package v3
