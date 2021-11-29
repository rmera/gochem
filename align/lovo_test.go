/*
 * lovo_test.go
 *
 * Copyright 2021 Raul Mera Adasme <rauldotmeraatusachdotcl>
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

package align

import (
	"fmt"
	"testing"

	chem "github.com/rmera/gochem"
)

func TestLovo(Te *testing.T) {
	path := "../test" //not a portable test, but I wanted a more "realistic" case. I'll change it for a more self contained test before
	molname := "/test_align.pdb"
	trajfilename := "/test_align.xtc"

	mol, err := chem.PDBFileRead(path+molname, false)
	if err != nil {
		fmt.Println("There was an error!", err.Error())
		Te.Error(err)
	}
	trajname := path + trajfilename
	o := DefaultOptions()
	o.NMostRigid = 10
	o.SetRigidPercent(90, 306)
	o.Skip = 0
	o.Cpus = 2
	o.WriteTraj = path + "/aligned.dcd"
	ret, err := LOVOnMostRigid(mol, mol.Coords[0], trajname, o)
	if err != nil {
		Te.Error(err)
	}
	//	fmt.Println("Most rigid residues:", ret.Nmols)
	fmt.Println(ret.PyMOLSel())
	fmt.Println(ret.String())
}
