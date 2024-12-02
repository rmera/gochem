/*
 * interfaces.go, part of gochem.
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

import v3 "github.com/rmera/gochem/v3"

/*The plan is equate PDBs XTCs and in the future DCDs. One needs to separate the molecule methods between actual molecule methods, that requires atoms and coordinates, []atom methods, and  DenseMatrix
 * methods. Then one must implements objects for Xtc trajs that are more or less equivalent to molecules and set up an interface so many analyses can be carried out exactly the same from
 * multiPDB or XTC or (eventually) DCD*/

// Traj is an interface for any trajectory object, including a Molecule Object
type Traj interface {

	//Is the trajectory ready to be read?
	Readable() bool

	//reads the next frame and returns it as DenseMatrix if keep==true, or discards it if false
	//it can also fill the (optional) box with the box vectors, it present in the frame.
	Next(output *v3.Matrix, box ...[]float64) error

	//Returns the number of atoms per frame
	Len() int
}

// ConcTraj is an interface for a trajectory that can be read concurrently.
type ConcTraj interface {

	//Is the trajectory ready to be read?
	Readable() bool

	/*NextConc takes a slice of bools and reads as many frames as elements the list has
	form the trajectory. The frames are discarted if the corresponding elemetn of the slice
	is false. The function returns a slice of channels through each of each of which
	a *matrix.DenseMatrix will be transmited*/
	NextConc(frames []*v3.Matrix) ([]chan *v3.Matrix, error)

	//Returns the number of atoms per frame
	Len() int
}

// Atomer is the basic interface for a topology.
type Atomer interface {

	//Atom returns the Atom corresponding to the index i
	//of the Atom slice in the Topology. Should panic if
	//out of range.
	Atom(i int) *Atom

	Len() int
}

// AtomChargerMultier is atomer but also gives a
// charge and multiplicity
type AtomMultiCharger interface {
	Atomer

	//Charge gets the total charge of the topology
	Charge() int

	//Multi returns the multiplicity of the topology
	Multi() int
}

// Masser can  return a slice with the masses of each atom in the reference.
type Masser interface {

	//Returns a column vector with the massess of all atoms
	Masses() ([]float64, error)
}

//Errors

//This error predates the "wrapping" error system of Go (i.e. the "%w" directive and the errors package). We should avoid
//using the Decorate method and/or make it use the "%w" directive internally.

// Error is the interface for errors that all packages in this library implement. The Decorate method allows to add and retrieve info from the
// error, without changing it's type or wrapping it around something else.
type Error interface {
	Error() string
	Decorate(string) []string //This is the new thing for errors. It allows you to add information when you pass it up. Each call also returns the "decoration" slice of strins resulting from the current call. If passed an empty string, it should just return the current value, not add the empty string to the slice.
	//The decorate slice should contain a list of functions in the calling stack, plus, for each function any relevant information, or nothing. If information is to be added to an element of the slice, it should be in this format: "FunctionName: Extra info"
}

// TrajError is the nterface for errors in trajectories
type TrajError interface {
	Error
	Critical() bool
	FileName() string
	Format() string
}

// LastFrameError has a useless function to distinguish the harmless errors (i.e. last frame) so  they can be
// filtered in a typeswith that looks for this interface.
type LastFrameError interface {
	TrajError
	NormalLastFrameTermination() //does nothing, just to separate this interface from other TrajError's

}
