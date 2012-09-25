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



import  "github.com/skelterjohn/go.matrix"


/*The plan is equate PDBs XTCs and in the future DCDs. One needs to separate the molecule methods between actual molecule methods, that requires atoms and coordinates, []atom methods, and  DenseMatrix
 * methods. Then one must implements objects for Xtc trajs that are more or less equivalent to molecules and set up an interface so many analyses can be carried out exactly the same from
 * multiPDB or XTC or (eventually) DCD*/

//Traj is an interface for any trajectory object, including a Molecule Object
type Traj interface{
	//Opens the file and prepares for reading, should also take care of the closing.
	Readable() bool
	
	//reads the next frame and returns it as DenseMatrix if keep==true, or discards it if false
	Next(keep bool) *matrix.DenseMatrix
	
	/*NextConc takes a slice of bools and reads as many frames as elements the list has
	form the trajectory. The frames are discarted if the corresponding elemetn of the slice
	is false. The function returns a slice of channels through each of each of which 
	a *matrix.DenseMatrix will be transmited*/
	NextConc(frames []bool)([]chan *matrix.DenseMatrix, error)
	
	//Selected, given a slice of ints, returns a matrix.DenseMatrix
	//containing the coordinates of the atoms with the corresponding index.
	SomeCoords(clist []int) (*matrix.DenseMatrix,error)
	
	//Returns the number of atoms per frame
	Len() int
	}


//Reference is an interface for any description of the type of atoms in a molecule,
//i.e. every characteristic of them, except for the coordinates and b-factors.
type Reference interface{

	//Charge gets the total charge of the topology
	 Charge()int
	
	//Unpaired gets the number of unpaired electrons in the topology
	 Unpaired()int
	
	//SetCharge sets the total charge of the topology to i
	 SetCharge(i int)
	
	//SetUnpaired sets the number of unpaired electrons in the topology to i
	 SetUnpaired(i int)
	
	//Atom returns the Atom corresponding to the index i
	//of the Atom slice in the Topology. Panics if 
	//out of range.
	 Atom(i int) (*Atom)
		
	//SetAtom sets the (i+1)th Atom of the topology to aM.
	//Panics if out of range
	 SetAtom(i int,at *Atom)
	
	//AddAtom appends an atom at the end of the topology
	 AddAtom(at *Atom)
	
	//SelectAtoms, given a list of ints,  returns an array of the atoms with the
	//corresponding position in the molecule
	//Changes to these atoms affect the original molecule.
	 SomeAtoms(atomlist []int) ([]*Atom, error)
	
	//Returns a column vector with the massess of all atoms
	//this will be changed to a tion that takes a Reference interface.
	 MassCol() (*matrix.DenseMatrix,error)
	
	//Returns the number of atoms in the reference
	 Len() int
	
	}




	
