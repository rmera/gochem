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
	"strings"
)

//Reduce uses the Reduce program (Word, et. al. (1999) J. Mol. Biol. 285,
//1735-1747. For more information see http://kinemage.biochem.duke.edu)
//To protonate a protein and flip residues. It writes the report from Reduce
//to a file called Reduce.err in the current dir. The Reduce executable must
//be a file called "reduce" and be in the PATH.
func Reduce(mol Atomer, coords *VecMatrix, build int, report *os.File) (*Molecule, error) {
	/*Unfortunately, I need to write each pdb to be protonated to disk. I just couldnt possibly make it work passing a PDB string
	  to reduce with a pipe. For some reason I only got a part of the PDB. It is something that should be fixed.*/
	pdbname := "gochemreducetmp.pdb"
	err := PDBWrite(pdbname, coords, mol, nil)

	if err != nil {
		return nil, err
	}
	flip := "-NOFLIP"
	if build == 1 {
		flip = "-FLIP"
	} else if build > 1 {
		flip = "-BUILD"
	}
	reduce := exec.Command("reduce", flip, pdbname) // , "-")
	//	inp, err := reduce.StdinPipe()
	//	if err != nil {
	//		return nil, err
	//	}
	out, err := reduce.StdoutPipe()
	if err != nil {
		return nil, err
	}
	out2, err := reduce.StderrPipe()
	if err != nil {
		return nil, err
	}
	if err := reduce.Start(); err != nil {
		return nil, err
	}
	/*	//Failed attempt to transmit the data directly to reduce using a pipe.
		binp := bufio.NewWriter(inp)
		chainprev:=mol.Atom(0).Chain
		outline:=""
		wcoord:=EmptyVecs()
		for i:=0;i<mol.Len();i++{
			wcoord.VecView(coords,i)
			outline,chainprev,err=writePDBLine(mol.Atom(i),wcoord,0.0,chainprev)
			_, err = binp.WriteString(outline)
			if err != nil {
				return nil, err
			}
			fmt.Println("tostdin", outline)
		}
		inp.Close()
	*/
	bufiopdb := bufio.NewReader(out)
	if report == nil {
		report, err := os.Create("Reduce.log")
		if err != nil {
			return nil, err
		}
		defer report.Close()
	}
	repio := bufio.NewWriter(report)
	dashes := 0
	for {
		s, err := bufiopdb.ReadString('\n')
		if err != nil {
			return nil, err
		}
		if strings.Contains(s, "USER  MOD ---------------------------------------------------") {
			dashes++
			if dashes > 1 {
				break
			}
		}
		repio.WriteString(s)
	}
	mol2, err := pdbBufIORead(bufiopdb, false)
	if err != nil {
		return nil, err
	}
	if _, err = repio.ReadFrom(out2); err != nil {
		return mol2, err
	}
	if err = reduce.Wait(); err != nil && !strings.Contains(err.Error(), "exit status 1") {
		return mol2, err
	}
	os.Remove(pdbname)
	return mol2, nil
}
