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

package chem

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"strings"

	v3 "github.com/rmera/gochem/v3"
	//	"runtime" /////////
	//	"fmt" ///////////
)

// Reduce uses the Reduce program (Word, et. al. (1999) J. Mol. Biol. 285,
// 1735-1747. For more information see http://kinemage.biochem.duke.edu)
// To protonate a protein and flip residues. It writes the report from Reduce
// to a file called Reduce.err in the current dir. The Reduce executable can be given,
// in "executable", otherwise it will be assumed that it is a file called "reduce" and it is in the PATH.
func Reduce(mol Atomer, coords *v3.Matrix, build int, report *os.File, executable string) (*Molecule, error) {
	//runtime.GOMAXPROCS(2) ////////////////////
	flip := "-NOFLIP"
	if build == 1 {
		flip = "-FLIP"
	} else if build > 1 {
		flip = "-BUILD"
	}
	if executable == "" {
		executable = "reduce"
	}
	reduce := exec.Command(executable, flip, "-") // , pdbname)
	inp, err := reduce.StdinPipe()
	if err != nil {
		return nil, err
	}
	out, err := reduce.StdoutPipe()
	if err != nil {
		return nil, CError{err.Error(), []string{"exec.StdoutPipe", "Reduce"}}
	}
	out2, err := reduce.StderrPipe()
	if err != nil {
		return nil, CError{err.Error(), []string{"exec.StderrPipe", "Reduce"}}

	}

	if err := reduce.Start(); err != nil {
		return nil, CError{err.Error(), []string{"exec.Start", "Reduce"}}
	}

	PDBWrite(inp, coords, mol, nil)
	inp.Close()

	bufiopdb := bufio.NewReader(out)
	if report == nil {
		report, err = os.Create("Reduce.log")
		if err != nil {
			return nil, CError{err.Error(), []string{"os.Create", "Reduce"}}

		}
		defer report.Close()
	}
	repio := bufio.NewWriter(report) //The reports of Reduce, which are written mostly to StdErr but also to StdOut will be compiled here.
	var stder_ready = make(chan error)
	/**
	In the following lines I am forced to use a gorutine to read the StdErr, as I need its contents, and reduce
	does not like them to be read serially. Apparently each of them get stuck at some point if you don't read the other,
	hence, I need to do it concurrently.
	**/
	go func() {
		defer repio.Flush()
		if _, err := repio.ReadFrom(out2); err != nil {
			out2.Close()
			stder_ready <- CError{err.Error(), []string{"bufio.ReadFrom", "Reduce"}}
		}
		out2.Close()
		stder_ready <- nil
	}()

	//I'm forced to use this buffer and later write to the report file to prevent a data race with the gorutine
	//reading StdErr.
	reportBuffer := make([]string, 0, 200)
	dashes := 0
	for {
		s, err := bufiopdb.ReadString('\n')
		if err != nil {
			return nil, fmt.Errorf("Reduce: %w", err)
		}
		reportBuffer = append(reportBuffer, s)

		if strings.Contains(s, "USER  MOD ---------------------------------------------------") {
			dashes++
			if dashes > 1 {
				break
			}
		}
	}

	mol2, err := pdbBufIORead(bufiopdb, false)
	if err != nil {
		return nil, fmt.Errorf("Reduce: %w", err)
	}
	out.Close()
	err = <-stder_ready
	//Here we don't need the buffering so we write directly to report (the file) instead of to repio.
	for _, v := range reportBuffer {
		_, err := report.Write([]byte(v))
		if err != nil {
			return mol2, err
		}
	}
	if err != nil {
		return mol2, err
	}

	return mol2, nil

}
