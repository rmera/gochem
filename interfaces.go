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
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package chem

import "github.com/rmera/gochem/v3"

/*The plan is equate PDBs XTCs and in the future DCDs. One needs to separate the molecule methods between actual molecule methods, that requires atoms and coordinates, []atom methods, and  DenseMatrix
 * methods. Then one must implements objects for Xtc trajs that are more or less equivalent to molecules and set up an interface so many analyses can be carried out exactly the same from
 * multiPDB or XTC or (eventually) DCD*/

//Traj is an interface for any trajectory object, including a Molecule Object
type Traj interface {

	//Is the trajectory ready to be read?
	Readable() bool

	//reads the next frame and returns it as DenseMatrix if keep==true, or discards it if false
	Next(output *v3.Matrix) error

	//Returns the number of atoms per frame
	Len() int
}

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

type Atomer interface {

	//Atom returns the Atom corresponding to the index i
	//of the Atom slice in the Topology. Should panic if
	//out of range.
	Atom(i int) *Atom

	Len() int
}

type ReadRef interface { /*NOTE: This method can probably be removed, as it doesn't seem to be used by any function that asks for these methods.*/
	Atomer

	//Returns a column vector with the massess of all atoms
	Masses() ([]float64, error)

	//Charge gets the total charge of the topology
	Charge() int

	//Multi returns the multiplicity of the topology
	Multi() int

	//Puts the atoms of Ref that are marked in the list of ints in the received.
	//Changes to these atoms affect the original molecule.
	//The charge and multiplicity (unpaired electrons) for the molecule is just the one
	//for the parent reference and its not guarranteed to be correct.
	SomeAtoms(Atomer, []int)
}

type WriteRef interface { /*NOTE: This whole interface can probably be removed, It appears that it is not used*/

	//Copy atoms in A into the received. This is a deep copy, so the received must have at least as many atoms as A
	CopyAtoms(A Atomer)

	//SetCharge sets the total charge of the topology to i
	SetCharge(i int)

	//SetUnpaired sets the number of unpaired electrons in the topology to i
	SetMulti(i int)

	//SetAtom sets the (i+1)th Atom of the topology to atm.
	//Panics if out of range
	SetAtom(i int, at *Atom)

	//AddAtom appends an atom at the end of the Ref
	AppendAtom(at *Atom)

	//Returns a copy of the Ref with the atom i deleted
	DelAtom(i int)

	//Changes the Ids and Molids of atoms for ones matching their current order
	ResetIDs()
}

//Ref (reference) is an interface for any description of the type of atoms in a molecule,
//i.e. every characteristic of them, except for the coordinates and b-factors.
//Read-write reference
type Ref interface { /*NOTE: WriteRef may be removed, which means that Ref would be removed too.*/
	ReadRef
	WriteRef
}

//Errors

type Error interface {
	Error() string
	Decorate(string) []string //This is the new thing for errors. It allows you to add information when you pass it up. Each call also returns the "decoration" slice of strins resulting from the current call. If passed an empty string, it should just return the current value, not add the empty string to the slice.
	//The decorate slice should contain a list of functions in the calling stack, plus, for each function any relevant information, or nothing. If information is to be added to an element of the slice, it should be in this format: "FunctionName: Extra info"
}

type TrajError interface {
	Error
	Critical() bool
	FileName() string
	Format() string
}

//Errors
type LastFrameError interface {
	TrajError
	NormalLastFrameTermination() //does nothing, just to separate this interface from other TrajError's

}
