/*
 * plot_test.go
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
 */

/*This provides some tests for the library functions requiring plotinum, in the form of little functions
 * that have practical applications*/

package chemplot

import (
	"fmt"
	"github.com/rmera/gochem"
	"testing"
)

//TestRama tests the Ramachandran plot functionality.
//it generates a Ramachandran plot for the chain A of the PDB 2c9v.
func TestRama(Te *testing.T) {
	mol, err := chem.PDBRead("../test/2c9v.pdb", true)
	if err != nil {
		Te.Error(err)
	}
	ramalist, err := RamaList(mol, "A", []int{0, -1}) ////
	if err != nil {
		Te.Error(err)
	}
	ramalist2, index := RamaResidueFilter(ramalist, []string{"HIS", "GLY"}, true)
	rama, err := RamaCalc(mol.Coords[0], ramalist2)
	if err != nil {
		Te.Error(err)
	}
	var i int
	for i = 0; i < len(ramalist); i++ {
		if index[i] != -1 {
			break
		}
	}
	fmt.Println("Rama", rama, len(rama), len(ramalist), mol.Len())
	err = RamaPlot(rama, []int{index[i]}, "Test Ramachandran", "../test/Rama")
	if err != nil {
		Te.Error(err)
	}
	//PDBWrite(mol,"test/Used4Rama.pdb")
	//for the 3 residue  I should get -131.99, 152.49.
}
