// +build puck

/*
 * pluckering_test.go
 *
 * Copyright 2013  <rmera@Holmes>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

package chem

//import "github.com/skelterjohn/go.matrix"
import "fmt"

//import "time"
import "testing"

//import "os"

//TestMultiXYZ tests that multi-XYZ files are opened and read correctly.
func TestXYZIO(Te *testing.T) {
	mol, err := XYZRead("test/5ring.xyz")
	if err != nil {
		fmt.Println("There was an error!")
		Te.Error(err)
	}
	newring:=Move5ring(Deg2Rad(40),0.5,mol.Coords[0])
	XYZWrite("test/new5ring2.xyz", mol, newring)
}
