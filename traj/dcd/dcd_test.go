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

package dcd

import (
	"fmt"
	"testing"

	chem "github.com/rmera/gochem"

	v3 "github.com/rmera/gochem/v3"
)

func TestDCDWrite2(Te *testing.T) {
	mol, traj, err := chem.XYZFileAsTraj("../../test/traj.xyz")
	if err != nil {
		fmt.Println("There was an error!", err.Error())
		Te.Error(err)
	}
	trajW, err := NewWriter("../../test/testW2.dcd", mol.Len())
	if err != nil {
		Te.Error(err)
	}

	for i := 0; ; i++ {
		fmt.Println("frame", i, "!!")
		err := traj.Next(mol.Coords[0])
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
		trajW.WNext(mol.Coords[0])
	}

	fmt.Println("XYZ read!")
}

//Tests the writing capabilities.
func TestDCDWrite(Te *testing.T) {
	fmt.Println("First test!")
	mol, err := chem.XYZFileRead("../../test/traj.xyz")
	if err != nil {
		Te.Error(err)
	}
	traj, err := NewWriter("../../test/testW2.dcd", mol.Len())
	if err != nil {
		Te.Error(err)
	}
	i := 0
	mat := v3.Zeros(mol.Len())
	for ; ; i++ {
		err := mol.Next(mat)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
		traj.WNext(mat)
	}
	fmt.Println("Over! frames read and written:", i)
}

/*TestDCD reads the frames of the test xtc file using the
 * "interactive" or "low level" functions, i.e. one frame at a time
 * It prints the firs 2 coordinates of each frame and the number of
 * read frames at the end.*/
func TestDCD(Te *testing.T) {
	fmt.Println("Fist test!")
	traj, err := New("../../test/test_align.dcd")
	if err != nil {
		Te.Error(err)
	}
	i := 0
	box := make([]float64, 9)
	mat := v3.Zeros(traj.Len())
	for ; ; i++ {
		err := traj.Next(mat, box)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			Te.Error(err)
			break
		}
		fmt.Println(mat.VecView(2), box)
		mat.Zero()
	}
	fmt.Println("Over! frames read:", i)
}

func TestFrameDCDConc(Te *testing.T) {
	traj, err := New("../../test/test_align.dcd")
	if err != nil {
		Te.Error(err)
	}
	frames := make([]*v3.Matrix, 3, 3)
	for i, _ := range frames {
		frames[i] = v3.Zeros(traj.Len())
	}
	results := make([][]chan *v3.Matrix, 0, 0)
	for i := 0; ; i++ {
		results = append(results, make([]chan *v3.Matrix, 0, len(frames)))
		coordchans, err := traj.NextConc(frames)
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				fmt.Printf("Frames read: %d\n", i*len(frames))
				break
			}
			Te.Error(err)
			break
		}
		for key, channel := range coordchans {
			results[len(results)-1] = append(results[len(results)-1], make(chan *v3.Matrix))
			go SecondRow(channel, results[len(results)-1][key], len(results)-1, key)
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

/*	for framebunch, j := range results {
		if j == nil {
			break
		}
		for frame, k := range j {
			if k == nil {
				fmt.Println(framebunch, frame)
				continue
			}
			fmt.Println(framebunch, frame, <-k)
		}
	}
}
*/
func SecondRow(channelin, channelout chan *v3.Matrix, current, other int) {
	if channelin != nil {
		temp := <-channelin
		viej := v3.Zeros(1)
		vector := temp.VecView(2)
		viej.Copy(vector)
		fmt.Println("sending througt", channelin, channelout, viej, current, other)
		channelout <- vector
	} else {
		channelout <- nil
	}
	return
}

/*
//TestFrameXtc reads the frames of the test xtc file from the first to
// the forth frame skipping one frame for each read one. It uses the
// "high level" function. It prints the frames read twince, and the
// coordinates of the forth atom of the last read frame
func TestFrameDCD(Te *testing.T) {
	fmt.Println("Second test!")
	traj,err:=NewDCD("test/test.dcd")
	if err != nil {
		Te.Error(err)
	}
	Coords, read, err := ManyFrames(traj,0, 5, 1)
	if err != nil {
		Te.Error(err)
	}
	fmt.Println(len(Coords), read, VecView(Coords[read-1],4))
	fmt.Println("DCD second test over!")
}

func TestFrameDCDConc(Te *testing.T) {
	fmt.Println("Third test!")
	traj,err:=NewDCD("test/test.dcd")
	if err != nil {
		Te.Error(err)
	}
	frames := []bool{true, true, true}
	results := make([][]chan *VecMatrix, 0, 0)
	_ = matrix.ZeroVecs(3, 3) //////////////
	for i := 0; ; i++ {
		results = append(results, make([]chan *VecMatrix, 0, len(frames)))
		coordchans, err := traj.NextConc(frames)
		if err != nil && err.Error() != "No more frames" {
			Te.Error(err)
		} else if err != nil {
			if coordchans == nil {
				break
			}
		}
		for key, channel := range coordchans {
			results[len(results)-1] = append(results[len(results)-1], make(chan *VecMatrix))
			go SecondRow(channel, results[len(results)-1][key], len(results)-1, key)
		}
	}
	for framebunch, j := range results {
		if j == nil {
			break
		}
		for frame, k := range j {
			if k == nil {
				fmt.Println(framebunch, frame)
				continue
			}
			fmt.Println(framebunch, frame, <-k)
		}
	}
}
*/
/*
func SecondRow(channelin, channelout chan *VecMatrix, current, other int) {
	if channelin != nil {
		temp := <-channelin
		vector := VecView(temp, 2)
		fmt.Println("sending througt", channelin, channelout, vector, current, other)
		channelout <- vector
	} else {
		channelout <- nil
	}
	return
}
*/
