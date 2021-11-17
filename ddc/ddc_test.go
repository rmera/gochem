/*
 * dcd_test.go
 *
 * Copyright 2012 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
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
 */
/*
 *
 *
 */

package ddc

import (
	"fmt"
	"log"
	"os"
	"testing"

	chem "github.com/rmera/gochem"

	v3 "github.com/rmera/gochem/v3"
)

/*TestDCD reads the frames of the test xtc file using the
 * "interactive" or "low level" functions, i.e. one frame at a time
 * It prints the firs 2 coordinates of each frame and the number of
 * read frames at the end.*/
func TestDDC(Te *testing.T) {
	fmt.Println("Fist test!")
	mol, err := chem.PDBFileRead("./test/weaita.pdb", true)
	if err != nil {
		Te.Error(err)
	}

	traj, err := New("./test/subset#000000")
	if err != nil {
		Te.Error(err)
	}
	fmt.Printf("Number of atoms in mol %d and traj %d \n", mol.Len(), traj.Len())
	if mol.Len() != traj.Len() {
		if err != nil {
			Te.Error(fmt.Errorf("Number of atoms in mol %d and traj %d should be the same\n", mol.Len(), traj.Len()))
		}

	}
	mat := v3.Zeros(traj.Len())
	err = traj.Next(mat)
	if err != nil {
		if _, ok := err.(chem.LastFrameError); ok {
			log.Print("EOF!")
		}
		Te.Error(err)
	}
	mol.Coords = append(mol.Coords, mat)
	out, err := os.Create("traj.pdb")
	if err != nil {
		Te.Error(err)
	}
	chem.MultiPDBWrite(out, mol.Coords, mol, nil)

}
