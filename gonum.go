/*
 * gonum.go, part of gochem.
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

//Package chem provides atom and molecule structures, facilities for reading and writing some
//files used in computational chemistry and some functions for geometric manipulations and shape
//indicators.
package chem

import "fmt"
import "github.com/skelterjohn/go.matrix"


/*Here I make a -very incomplete- implementation of the gonum api backed by go.matrix, which will enable me to port gochem to gonum. 
 * Since the agreement in the gonum community was NOT to build a temporary implementation, I just build the functions that
 * gochem uses, only the type Float64, and do not export any of the functions.
 * all the names will start with gn (i.e. RandomFunc becomes gnRandomFunc) so its latter easy to use search and replace to set the 
 * correct import path when gonum is implemented (such as gonum.RandomFunc)*/
 
 
 
type gnFloat64 matrix.DenseMatrix
 
 
 func gnNewFloat64(data []float64,rows,cols int) *gnFloat64{
	return gnFloat64(matrix.MakeDenseMatrix(data, rows, cols))
}
  

(F *gnFloat64) Dims()(int, int){
	return F.Rows(),F.Cols()
}

 
 
 
 
 
 
 
 
 
 
 
 
 
 
