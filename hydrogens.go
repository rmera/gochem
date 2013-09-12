/*
 * hydrogens.go, part of gochem.
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

import (
	"bufio"
	"os"
	"os/exec"
)

//Reduce uses the Reduce program
//(Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.
//For more information see http://kinemage.biochem.duke.edu)
//To protonate a protein and flip residues. It writes
//the report from Reduce to a file called Reduce.err in the current dir.
//The Reduce executable must be a file called "reduce" and be in the PATH.
func Reduce(mol Atomer, coords *CoordMatrix, build int, report string) (*Molecule, error) {
	pdb, err := PDBStringWrite(mol, coords, nil)
	if err != nil {
		return nil, err
	}
	flip := "-NOFLIP"
	if build == 1 {
		flip = "-FLIP"
	} else if build > 1 {
		flip = "-BUILD"
	}
	reduce := exec.Command("reduce", flip, "-")
	inp, err := reduce.StdinPipe()
	if err != nil {
		return nil, err
	}
	out, err := reduce.StdoutPipe()
	if err != nil {
		return nil, err
	}
	out2, err := reduce.StderrPipe()
	if err != nil {
		return nil, err
	}
	defer out.Close()
	defer out2.Close()
	if err := reduce.Start(); err != nil {
		return nil, err
	}
	binp := bufio.NewWriter(inp)
	_, err = binp.WriteString(pdb)
	if err != nil {
		return nil, err
	}
	inp.Close()
	bufiopdb := bufio.NewReader(out)
	mol2, err := pdbBufIORead(bufiopdb, false)
	rep, err := os.Create(report)
	if err != nil {
		return mol2, err
	}
	defer rep.Close()
	repio := bufio.NewWriter(rep)
	_, err = repio.ReadFrom(out2)
	if err != nil {
		return mol2, err
	}
	return mol2, nil
}
