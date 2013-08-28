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

/*The plan is equate PDBs XTCs and in the future DCDs. One needs to separate the molecule methods between actual molecule methods, that requires atoms and coordinates, []atom methods, and  DenseMatrix
 * methods. Then one must implements objects for Xtc trajs that are more or less equivalent to molecules and set up an interface so many analyses can be carried out exactly the same from
 * multiPDB or XTC or (eventually) DCD*/

//Traj is an interface for any trajectory object, including a Molecule Object
type Traj interface {

	//Is the trajectory ready to be read?
	Readable() bool

	//reads the next frame and returns it as DenseMatrix if keep==true, or discards it if false
	Next(output *CoordMatrix) error

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
	NextConc(frames []*CoordMatrix) ([]chan *CoordMatrix, error)

	//Returns the number of atoms per frame
	Len() int
}

//Ref (reference) is an interface for any description of the type of atoms in a molecule,
//i.e. every characteristic of them, except for the coordinates and b-factors.
type Ref interface {

	//Charge gets the total charge of the topology
	Charge() int

	//Unpaired gets the number of unpaired electrons in the topology
	Unpaired() int

	//SetCharge sets the total charge of the topology to i
	SetCharge(i int)

	//SetUnpaired sets the number of unpaired electrons in the topology to i
	SetUnpaired(i int)

	//Atom returns the Atom corresponding to the index i
	//of the Atom slice in the Topology. Should panic if
	//out of range.
	Atom(i int) *Atom

	//SetAtom sets the (i+1)th Atom of the topology to atm.
	//Panics if out of range
	SetAtom(i int, at *Atom)

	//AddAtom appends an atom at the end of the Ref
	AddAtom(at *Atom)

	//Puts the atoms of Ref that are marked in the list of ints in the received.
	//Changes to these atoms affect the original molecule.
	//The charge and multiplicity (unpaired electrons) for the molecule is just the one
	//for the parent reference and its not guarranteed to be correct.
	SomeAtoms(Ref, []int)

	//Returns a column vector with the massess of all atoms
	//this will be changed to a tion that takes a Reference interface.
	MassCol() (*CoordMatrix, error)

	//Returns a copy of the Ref with the atom i deleted
	DelAtom(i int)

	//Changes the Ids and Molids of atoms for ones matching their current order
	ResetIds()

	//Returns the number of atoms in the reference
	Len() int
}

/*
//This allows to set QM calculations using different programs.
//Currently ORCA and MOPAC (2009/2012) are supported.
type QMRunner interface {
	//Sets the number of  CPUs for the calculation, when possible
	SetnCPU(cpu int)

	//Sets the name for the job, used for input
	//and output files. The extentions will depend on the program.
	SetName(name string)

	//Sets the command to run the QM  program.
	SetCommand(name string)

	//Sets the defaults for a QM calculation. Varies depending on the
	// program,
	SetDefaults()

	//BuildInput builds an input for the QM program based int the data in
	//atoms, coords and C. returns only error.
	BuildInput(atoms Ref, coords *matrix.DenseMatrix, Q *QMCalc) error

	//Run runs the QM program for a calculation previously set.
	//it waits or not for the result depending of the value of
	//wait.
	Run(wait bool) (err error)

	//GetEnergy gets the last energy for a  calculation by parsing the
	//QM program's output file. Return error if fail. Also returns
	//Error ("Probable problem in calculation")
	//if there is a energy but the calculation didnt end properly.
	GetEnergy() (float64, error)

	//Get Geometry reads the optimized geometry from a calculation
	//output. Returns error if fail. Returns Error ("Probable problem
	//in calculation") if there is a geometry but the calculation didnt
	//end properly*
	GetGeometry(atoms Ref) (*matrix.DenseMatrix, error)
}
*/
