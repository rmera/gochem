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

package csf

import (
	"fmt"
	"testing"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/traj/dcd"
	"github.com/rmera/gochem/traj/xtc"
	v3 "github.com/rmera/gochem/v3"
)

//Tests the writing capabilities.
func TestWrite(Te *testing.T) {
	fmt.Println("STF write test!")
	mol, err := chem.PDBFileRead("../../test/test_stf.pdb", false)
	if err != nil {
		Te.Error(err)
	}
	rtraj, err := xtc.New("../../test/test_stf.xtc")
	if err != nil {
		Te.Error(err)
	}
	wtraj, err := NewWriter("../../test/test_stf.ctz", mol, nil, 9)
	if err != nil {
		Te.Error(err)
	}
	defer wtraj.Close()
	i := 0
	mat := v3.Zeros(mol.Len())
	for ; ; i++ {
		err := rtraj.Next(mat)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
		wtraj.WNext(mat)
	}
	fmt.Println("Over! frames read and written:", i)
}

//Now the read
func TestRead(Te *testing.T) {
	fmt.Println("STF read test!")

	mol, err := chem.PDBFileRead("../../test/test_stf.pdb", false)
	if err != nil {
		Te.Error(err)
	}
	rtraj, _, err := New("../../test/test_stf.ctz")
	if err != nil {
		Te.Error(err)
	}
	wtraj, err := dcd.NewWriter("../../test/test_stf.dcd", mol.Len())
	if err != nil {
		Te.Error(err)
	}
	defer rtraj.Close()
	i := 0
	mat := v3.Zeros(rtraj.Len())
	for ; ; i++ {
		err := rtraj.Next(mat)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
		wtraj.WNext(mat)
	}
	fmt.Println("Over! frames read:", i)
}

func BenchmarkWriteDCD(B *testing.B) {
	fmt.Println("DCD write bench!")
	mol, err := chem.PDBFileRead("../../test/test_stf.pdb", false)
	if err != nil {
		B.Error(err)
	}
	rtraj, err := xtc.New("../../test/test_stf.xtc")
	if err != nil {
		B.Error(err)
	}
	wtraj, err := dcd.NewWriter("../../test/test_b.dcd", mol.Len())
	if err != nil {
		B.Error(err)
	}
	i := 0
	mat := v3.Zeros(mol.Len())
	for ; ; i++ {
		err := rtraj.Next(mat)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			B.Error(err)
			break
		}
		wtraj.WNext(mat)
	}
}

func BenchmarkWriteSTF(B *testing.B) {
	fmt.Println("STF write bench!")
	mol, err := chem.PDBFileRead("../../test/test_stf.pdb", false)
	if err != nil {
		B.Error(err)
	}
	rtraj, err := xtc.New("../../test/test_stf.xtc")
	if err != nil {
		B.Error(err)
	}
	wtraj, err := NewWriter("../../test/test_b.stf", mol, nil)
	if err != nil {
		B.Error(err)
	}
	i := 0
	mat := v3.Zeros(mol.Len())
	for ; ; i++ {
		err := rtraj.Next(mat)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			B.Error(err)
			break
		}
		wtraj.WNext(mat)
	}
}

//Now the read
func BenchmarkReadSTF(Te *testing.B) {
	fmt.Println("STF read test!")

	//	mol, err := chem.PDBFileRead("test.pdb", false)
	//	if err != nil {
	//		Te.Error(err)
	//	}
	rtraj, _, err := New("../../test/test.ctz")
	if err != nil {
		Te.Error(err)
	}
	defer rtraj.Close()
	i := 0
	mat := v3.Zeros(rtraj.Len())
	for ; ; i++ {
		err := rtraj.Next(mat)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
	}
	fmt.Println("Over! frames read:", i)
}
