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

package stf

import (
	"fmt"

	"testing"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/dcd"
	"github.com/rmera/gochem/xtc"

	v3 "github.com/rmera/gochem/v3"
)

//Tests the writing capabilities.
func TestSTFWrite(Te *testing.T) {
	var err error
	fmt.Println("STF write test!")
	mol, err := chem.PDBFileRead("../test/test_stf.pdb", false)
	if err != nil {
		Te.Error(err)
	}
	rtraj, err := xtc.New("../test/test_stf.xtc")
	if err != nil {
		Te.Error(err)
	}
	wtraj, err := NewWriter("../test/test_stf.stz", mol.Len(), nil)
	if err != nil {
		Te.Error(err)
	}
	defer wtraj.Close()
	i := 0
	mat := v3.Zeros(mol.Len())
	for ; ; i++ {
		if i%101 == 0 {
			err = rtraj.Next(mat)
		} else {
			err = rtraj.Next(nil)
		}
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
		if i%101 == 0 {
			wtraj.WNext(mat)
		}
	}
	fmt.Println("Over! frames read and written:", i)
}

//Now the read
func TestSTF(Te *testing.T) {
	fmt.Println("STF read test!")

	mol, err := chem.PDBFileRead("../test/test_stf.pdb", false)
	if err != nil {
		Te.Error(err)
	}
	rtraj, _, err := New("../test/test_stf.stz")
	if err != nil {
		Te.Error(err)
	}
	wtraj, err := dcd.NewWriter("../test/test_stf.dcd", mol.Len())
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

func TestConc(Te *testing.T) {
	fmt.Println("Concurrency test!")
	traj, _, err := New("../test/test_stf.stz")
	if err != nil {
		Te.Error(err)
	}
	frames := make([]*v3.Matrix, 3, 3)
	for i, _ := range frames {
		frames[i] = v3.Zeros(traj.Len())
	}
	frames[1] = nil /////Just a test
	results := make([][]chan *v3.Matrix, 0, 0)
	for i := 0; ; i++ {
		results = append(results, make([]chan *v3.Matrix, 0, len(frames)))
		coordchans, err := traj.NextConc(frames)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
		for key, channel := range coordchans {
			results[len(results)-1] = append(results[len(results)-1], make(chan *v3.Matrix))
			go LastRow(channel, results[len(results)-1][key], len(results)-1, key)
		}
		res := len(results) - 1
		for frame, k := range results[res] {
			if k == nil {
				fmt.Println(frame)
				continue
			}
			fmt.Println(res, frame, <-k)
		}
	}
}

func LastRow(channelin, channelout chan *v3.Matrix, current, other int) {
	var vector *v3.Matrix
	if channelin != nil {
		temp := <-channelin
		viej := v3.Zeros(1)
		if temp != nil {
			vector = temp.VecView(temp.Len() - 1)
			viej.Copy(vector)
		}
		fmt.Println("sending througt", channelin, channelout, viej, current, other)
		channelout <- vector
	} else {
		channelout <- nil
	}
	return
}

func BenchmarkWriteDCD(B *testing.B) {
	fmt.Println("DCD write bench!")
	mol, err := chem.PDBFileRead("../test/test_stf.pdb", false)
	if err != nil {
		B.Error(err)
	}
	rtraj, err := xtc.New("../test/test_stf.xtc")
	if err != nil {
		B.Error(err)
	}
	wtraj, err := dcd.NewWriter("../test/test_stf_b.dcd", mol.Len())
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
	mol, err := chem.PDBFileRead("../test/test_stf.pdb", false)
	if err != nil {
		B.Error(err)
	}
	rtraj, err := xtc.New("../test/test_stf.xtc")
	if err != nil {
		B.Error(err)
	}
	wtraj, err := NewWriter("../test/test_stf_b.stz", mol.Len(), nil)
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
