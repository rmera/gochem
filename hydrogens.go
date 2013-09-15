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
	"fmt"
//	"unicode/utf8"
//	"strings"
)

//Reduce uses the Reduce program
//(Word, et. al. (1999) J. Mol. Biol. 285, 1735-1747.
//For more information see http://kinemage.biochem.duke.edu)
//To protonate a protein and flip residues. It writes
//the report from Reduce to a file called Reduce.err in the current dir.
//The Reduce executable must be a file called "reduce" and be in the PATH.
func Reduce(mol Atomer, coords *CoordMatrix, build int, report *bufio.Writer) (*Molecule, error) {
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
	if err := reduce.Start(); err != nil {
		return nil, err
	}
	binp := bufio.NewWriter(inp)
	l:=len([]byte(pdb))
	l2, err := binp.WriteString(pdb)
	if err != nil {
		return nil, err
	}
	fmt.Println(pdb, "the full pdb")
	fmt.Println(l,l2)
	inp.Close()
	if report==nil{
		rep, err := os.Create("Reduce.log")
		if err != nil {
			return nil, err
		}
		defer rep.Close()
		report = bufio.NewWriter(rep)
	}
	//the first lines of the PDB belong to the report.
/*	for {
		fmt.Println("writing lines!")
		line,err:=bufiopdb.ReadString('\n')
		if err!=nil{
			return nil, err
		}
		if strings.Contains(line, "USER  MOD ----------------------------------------------------------"){
			break
		}
		report.WriteString(line) //NO error checking!!!
	}
*/
	fmt.Println("will write")
	mol2, err := PDBReaderRead(out, false)
	_, err = report.ReadFrom(out2)
	if err != nil {
		return mol2, err
	}
	fmt.Println("wrote")
	if err := reduce.Wait(); err != nil && err.Error()!="exit status 1"{
		return nil, err
	}
	fmt.Println("last", mol2.Atom(mol2.Len()-1))
	return mol2, nil
}
