// +build dcd

/*
 * untitled.go
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package chem

import "fmt"
import "testing"

/*TestXtc reads the frames of the test xtc file using the
 * "interactive" or "low level" functions, i.e. one frame at a time
 * It prints the firs 2 coordinates of each frame and the number of
 * read frames at the end.*/
func estDCD(Te *testing.T) {
	fmt.Println("Fist test!")
	traj, err := NewDCD("test/test.dcd")
	if err != nil {
		Te.Error(err)
	}
	i := 0
	mat := ZeroVecs(traj.Len())
	for ; ; i++ {
		err := traj.Next(mat)
		if err != nil && err.Error() != "No more frames" {
			Te.Error(err)
			break
		} else if err == nil {
			fmt.Println(VecView(mat, 2))
		} else {
			break
		}
	}
	fmt.Println("Over! frames read:", i)
}

func TestFrameDCDConc(Te *testing.T) {
	traj, err := NewDCD("test/test.dcd")
	if err != nil {
		Te.Error(err)
	}
	frames := make([]*VecMatrix, 3, 3)
	for i, _ := range frames {
		frames[i] = ZeroVecs(traj.Len())
	}
	results := make([][]chan *VecMatrix, 0, 0)
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
func SecondRow(channelin, channelout chan *VecMatrix, current, other int) {
	if channelin != nil {
		temp := <-channelin
		viej := ZeroVecs(1)
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
