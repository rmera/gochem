/*
 * basic.go, part of gochem.
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/


//Package chem provides atom and molecule structures, facilities for reading and writing some
//files used in computational chemistry and some functions for geometric manipulations and shape
//indicators.
package chem

import "fmt"
import  "github.com/skelterjohn/go.matrix"









//Atom contains the atoms read except for the coordinates, which will be in a matrix
//and the b-factors, which are in a separate slice of float64.
type Atom struct{
	Name string
	Id int
	Tag int //Just added this for something that someone might want to keep that is not a float.
	Molname string
	Molname1 byte   //the one letter name for residues and nucleotids
	Molid int
	Chain byte
	Mass float64 //hopefully all these float64 are not too much memory
	Occupancy float64
	Vdw float64 
	Charge float64
	Symbol string
	Het bool  // is hetatm in the pdb file?
	}

//Molecule contains all the info for a molecule in many states. The info that is expected to change between states,
//Coordinates and b-factors are stored separately from other atomic info.
type Molecule struct{
	Atoms []*Atom
	Coords []*matrix.DenseMatrix
	Bfactors [][]float64
	Total_charge int
	Unpaired_electrons int
	}

//The molecule methods:

//GetMassArray return an array with massess of atoms and an error if they have not been calculated
func (M *Molecule) GetMassArray() ([]float64,error){
	if M==nil{
		return nil,fmt.Errorf("Molecule is nil") 
		}
	if M.Atoms==nil{
		return nil, fmt.Errorf("Molecule is empty")
		}
	mass:=make([]float64,len(M.Atoms))
	for i := range M.Atoms{
		if M.Atoms[i].Mass==0{
			return nil, fmt.Errorf("Not all the masses have been obtained")
			}
		mass[i]=M.Atoms[i].Mass
		}
	
	return mass, nil
	}
	
//AddFrame akes a matrix of coordinates and appends them at the end of the Coords.
// It checks that the number of coordinates matches the number of atoms.
func (M *Molecule) AddFrame(newframe *matrix.DenseMatrix) error{
	if M==nil{
		return fmt.Errorf("Molecule is nil") 
		}
	if M.Atoms==nil{
		return fmt.Errorf("Molecule is empty")
		}
	if newframe == nil{
		return fmt.Errorf("Attempted to add nil frame") 
		}
	if newframe.Cols()!=3{
		return fmt.Errorf("Malformed coord matrix!") 
		}
	if len(M.Atoms)!=newframe.Rows(){
		//fmt.Println("lens:", len(M.Atoms),len(newframe)/3)
		return fmt.Errorf("Wrong number of coordinates in frame!")  //change this for something that prints the number of atoms and coords
		}
	if M.Coords==nil{
		M.Coords=make([]*matrix.DenseMatrix,1,1)
		}
	M.Coords=append(M.Coords,newframe)
	return nil
	}

//AddManyFrames adds the array of matrices newfames to the molecule. It checks that
//the number of coordinates matches the number of atoms.
func (M *Molecule) AddManyFrames(newframes []*matrix.DenseMatrix) error{
	if M==nil{
		return fmt.Errorf("Molecule is nil") 
		}
	if M.Atoms==nil{
		return fmt.Errorf("Molecule is empty")
		}
	if newframes == nil {
		return fmt.Errorf("Attempted to add nil frames") 
		}
	if M.Atoms == nil{
		return fmt.Errorf("The molecule has no atoms") 
		}
	atomslen:=len(M.Atoms)
	if M.Coords==nil{
		M.Coords=make([]*matrix.DenseMatrix,1,len(newframes))
		}
	//Must add something here to change
	for i:=range newframes{
		if atomslen!=newframes[i].Rows(){
			return fmt.Errorf("Wrong number of coordinates (%d)  in frame %d!",i,newframes[i].Rows()) 
			}
		M.Coords=append(M.Coords,newframes[i]) //At the end there shouldnt be so much copying since
	   //only pointers are copied.
		}
	return nil
	}


/*
//Gets you slice of coordinates from atom first to atom last in frame frame

func (M *Molecule) Coords(first int, last int, frame int) []float64{
	return M.coords[frame][first*3:(last*3)+3]
	}
*/


//Coord returns a DenseMatrix with the coordinates of the atom atom in the frame frame
//Changes to this DenseMatrix affect the original coordinates.
//Notice that this doesnt return errors at all, just nil if there is something wrong.
//This is because this function is meant to be embebed in bigger expressions. Check before using!
func (M *Molecule) Coord(atom, frame int) *matrix.DenseMatrix{
	if frame>=len(M.Coords) || atom >= M.Coords[0].Rows(){
		return nil
		}
	return M.Coords[frame].GetRowVector(atom)
	}

//GetCoordsSlice, given a list of ints, and a frame, returns an slice of DenseMatrix where the nth
//element contains the coordinates to the atom in the position clist[n] in M.Atoms.
//Changes to these matrices affect the original M.Coords. It checks for correctness of the frame and the
//Atoms requested.
func (M *Molecule) GetCoordsSlice(clist []int, frame int) ([]*matrix.DenseMatrix,error){
	if M==nil{
		return nil,fmt.Errorf("Molecule is nil") 
		}
	if M.Atoms==nil{
		return nil, fmt.Errorf("Molecule is empty")
		}
	var err error
	var ret []*matrix.DenseMatrix
	if frame>=len(M.Coords){
		return ret,fmt.Errorf("Frame requested (%d) out of range!",frame)
		}
	lencoords:=M.Coords[frame].Rows()
	for k,j:=range(clist){
		if j>lencoords-1{
			return ret,fmt.Errorf("Coordinate requested (Number: %d, value: %d) out of range!",k,j)
			}
		ret=append(ret,M.Coords[frame].GetRowVector(j))
		}
	return ret,err
	}
	
	
//GetCoords, given a list of ints and the desired frame, returns an slice matrix.DenseMatrix
//containing the coordinates of the atoms with the corresponding index.
//This function returns a copy, not a reference, so changes to the returned matrix
//don't alter the original. It check for correctness of the frame and the
//Atoms requested.
func (M *Molecule) GetCoords(clist []int, frame int) (*matrix.DenseMatrix,error){
	var err error
	var ret [][]float64
	if frame>=len(M.Coords){
		return nil,fmt.Errorf("Frame requested (%d) out of range!",frame)
		}
	lencoords:=M.Coords[frame].Rows()
	for k,j:=range(clist){
		if j>lencoords-1{
			return nil,fmt.Errorf("Coordinate requested (Number: %d, value: %d) out of range!",k,j)
			}
		tmp:=M.Coords[frame].GetRowVector(j).Array()
		if len(tmp)!=3{
			return nil,fmt.Errorf("Coordinate %d has %d components instead of 3",k, len(tmp))
			}
		ret=append(ret,tmp)
		}
	return matrix.MakeDenseMatrixStacked(ret),err
/*
 * OLD IMPLEMENTATION, STILL NOT SURE OF WHICH IS BETTER:
 * 	//I have to treat the first element separately in order to start with a filled ret to stack
	if clist[0]>lencoords-1{
		return nil,fmt.Errorf("Coordinate requested (Number: 0, value: %d) out of range!",clist[0])
		}
	ret:=M.Coords[frame].GetRowVector(clist[0])
	for k,j:=range(clist[1:]){
		if j>lencoords-1{
			return ret,fmt.Errorf("Coordinate requested (Number: %d, value: %d) out of range!",k,j) //This could say which element gave the problem
			}
		ret,err=ret.Stack(M.Coords[frame].GetRowVector(j))
		}
	return ret,err
	*/
	}


//GetAtoms, given a list of ints,  returns an array of the atoms with the
//corresponding position in the molecule
//Changes to these atoms affect the original molecule.
func (M *Molecule) GetAtoms(atomlist []int) ([]*Atom, error){
	var err error
	var ret []*Atom
	lenatoms:=len(M.Atoms)
	for k,j:=range(atomlist){
		if j>lenatoms-1{
			return nil,fmt.Errorf("Atom requested (Number: %d, value: %d) out of range!",k,j)
			}
		ret=append(ret,M.Atoms[j])
		}
	return ret,err
	}

//SetCoords replaces the coordinates of atoms in the positions given by atomlist with the ones in newcoords (in order)
//If atomlist contains a single element, it replaces as many coordinates as given in newcoords, starting 
//at the element in atomlist. In the latter case, the function checks that there are enough coordinates to
//replace and returns an error if not.
func (M *Molecule) SetCoords(atomlist []int, frame int, newcoords *matrix.DenseMatrix) (error){
	if frame>=len(M.Coords){
		return fmt.Errorf("Frame (%d) out of range!",frame)
		}
	//If supplies a list with one number, the newcoords will replace the old coords
	//Starting that number. We do check that you don't put more coords than spaces we have.
	if len(atomlist)==1{
		if newcoords.Rows()>M.Coords[frame].Rows()-atomlist[0]-1{
			return fmt.Errorf("Cant replace starting from position %d: Not enough atoms in molecule", atomlist[0])
			} 
		M.Coords[frame].SetMatrix(atomlist[0],0,newcoords)
		return nil
		}
	//If the list has more than one atom
	lenatoms:=len(M.Atoms)	
	for k,j:=range(atomlist){
		if j>lenatoms-1{
			return fmt.Errorf("Requested position number: %d (%d) out of range",k,j)
			}
		M.Coords[frame].SetMatrix(j,0,newcoords.GetRowVector(k))
		}
	return nil
	}
	
//Swap function, as demanded by sort.Interface. It swaps atoms, coordinates 
//(all frames) and bfactors of the molecule.	
func (M *Molecule) Swap(i,j int) {
	M.Atoms[i],M.Atoms[j]=M.Atoms[j],M.Atoms[i]
	for k:=0;k<len(M.Coords);k++{
		M.Coords[k].SwapRows(i,j)
		M.Bfactors[k][i],M.Bfactors[k][j]=M.Bfactors[k][j],M.Bfactors[k][i]
		}
	}

//Less: Should the atom i be sorted before atom j?
func (M *Molecule) Less (i, j int) bool {
	return M.Bfactors[0][i]<M.Bfactors[0][j]
	}

//Len return the length of the molecule.
func (M *Molecule) Len() int{
	return len(M.Atoms)
	}
	
//LenFrames returns the number of frames in the molecule
func (M *Molecule) LenFrames() int{
	return len(M.Coords)
	}
	
//Corrupted checks whether the molecule is corrupted, i.e. the
//coordinates don't match the number of atoms. It also checks
//That the coordinate matrices have 3 columns.
func (M *Molecule) Corrupted() error{
	var err error
	lastbfac:=len(M.Bfactors)-1
	for i:=range M.Coords{
		if len(M.Atoms)!=M.Coords[i].Rows() || M.Coords[i].Cols()!=3{
			//Not tested.
			err=fmt.Errorf("Inconsistent coordinates/atoms in frame %d", i) 
			break
			}
		//Since bfactors are not as important as coordinates, we will just fill with 
		//zeroes anything that is lacking or incomplete instead of returning an error.
		if lastbfac<i{
			bfacs:=make([]float64,len(M.Atoms),len(M.Atoms))
			M.Bfactors=append(M.Bfactors,bfacs)
			}else if len(M.Bfactors[i])<len(M.Atoms){
			bfacs:=make([]float64,len(M.Atoms),len(M.Atoms))
			M.Bfactors[i]=bfacs
			}
		}
	return err
	}
//End Molecule methods













