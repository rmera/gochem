// +build ignore

/*
 * init_goblas.go, part of gochem.
 *
 * Copyright 2014 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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

//The point of this is that the user can use compile tags to get their programs compiled with goblas or cblas.
//goblas is the default, the cblas tag being required to use, well, cblas.

package v3

import (
	"gonum.org/v1/gonum/blas/native"
	"gonum.org/v1/gonum/mat"
)

func init() {
	mat.Register(native.Blas{})
}
