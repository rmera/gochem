/*
 * Rama.go, part of gochem
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 * 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2.1 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 * 
 */
 /***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package chem
 
import (
		"fmt"
		"strings"
		"math"
        "code.google.com/p/plotinum/plot"
        "code.google.com/p/plotinum/plotter"
        "image/color"
)

func RamaPlot(fulldata [][][]float64, plotname string) error{
	for number,data:=range(fulldata){
		if data==nil{
			return fmt.Errorf("Given nil data")
			}
			pts := make(plotter.XYs, len(data)) //just first frame for now
			//this might not be too efficient
		for key,val:=range(data){ 
			pts[key].X=val[0]
			pts[key].Y=val[1]
			}
		// Create a new plot, set its title and
		// axis labels.
		p, err := plot.New()
			if err != nil {
				return err
				}
		p.Title.Text = "Rama plot"
		p.X.Label.Text = "Phi"
		p.Y.Label.Text = "Psi"
		//Constant axes
		p.X.Min=-180
		p.X.Max=180
		p.Y.Min=-180
		p.Y.Max=180
		// Draw the grid 
		p.Add(plotter.NewGrid())
		// Make a scatter plotter and set its style.
		s := plotter.NewScatter(pts)
		s.GlyphStyle.Color = color.RGBA{R: 255, A: 255}
		// Add the plotter 
		p.Add(s)
		// Save the plot to a PNG file.
		filename:=fmt.Sprintf("%s%d.png",plotname,number)
		if err := p.Save(4, 4, filename); err != nil {
			return err
		}
	}
	return nil
}








func RamaCalc(M *Molecule, dihedrals [][]int, frames []int)([][][]float64,error){
	if M == nil || dihedrals == nil || frames == nil{
		return nil,  fmt.Errorf("Given nil data")
		}
	if err:=M.Corrupted(); err!=nil{
		return nil, fmt.Errorf("Molecule corrupted") 
		}
	if len(M.Coords)<=frames[len(frames)-1]{
		return nil,  fmt.Errorf("Frame out of range") 
		} 	
	Rama:=make([][][]float64,0,len(frames))
	for _,i:=range(frames){
		Ramaframe:=make([][]float64,0,len(dihedrals))
		for id,j:=range(dihedrals){
			Cprev:=M.Coords[i].GetRowVector(j[0])
			N:=M.Coords[i].GetRowVector(j[1])
			Ca:=M.Coords[i].GetRowVector(j[2])
			C:=M.Coords[i].GetRowVector(j[3])
			Npost:=M.Coords[i].GetRowVector(j[4])
			phi,err:=Dihedral(Cprev,N,Ca,C)
			psi,err1:=Dihedral(N,Ca,C,Npost)
			if err!=nil || err1!=nil{
				return nil, fmt.Errorf("Dihedral calculation failed for residue %d",id+2) 
				}
			temp:=[]float64{phi*(180/math.Pi),psi*(180/math.Pi)}
			Ramaframe=append(Ramaframe,temp)
			}
		Rama=append(Rama,Ramaframe)
		}
	return Rama, nil
	}

//isIn is a helper for the RamaList function, 
//returns true if test is in container, false otherwise.
func IsInInt(container []int, test int) bool{
	if container==nil{
		return false
		}
	for _,i:=range(container){
		if test==i{
			return true
			}
		}
	return false
	}

/*RamaList takes a molecule and obtains a list of lists of five int. Each element
 * contain the indexes needed for one dihedral of a Rama plot. It gets the dihedral
 * indices for all atoms in the range resran, if resran has 2 elements defining the 
 * boundaries. Otherwise, returns dihedral lists for the residues included in 
 * resran. If resran has 2 elements and the last is -1, RamaList will
 * get all the dihedral for residues from resran[0] to the end of the chain.
 * It only obtain dihedral lists for residues belonging to a chain included in chains */
func RamaList(M *Molecule, chains string,resran []int) ([][]int, error){
	RamaList:=make([][]int,0,0)
	if len(resran)==2{
		if resran[1]==-1{
			resran[1]=999999999 //should work!
			}
	}
	if M == nil{
		return nil, fmt.Errorf("nil Molecule")
		}
	if err:=M.Corrupted(); err!=nil{
		return nil, fmt.Errorf("Molecule corrupted") //Must change this
		} 
	C:=-1
	N:=-1
	Ca:=-1
	Cprev:=-1
	Npost:=-1
	chainprev:=byte('9') //any non-valid chain name
	for num,at:=range(M.Atoms){
		//First get the indexes we need
		if strings.Contains(chains,string(at.Chain)) || at.Chain==' '{
			if at.Chain!=chainprev  {
				chainprev=at.Chain
				C=-1
				N=-1
				Ca=-1
				Cprev=-1
				Npost=-1
				}
			if at.Name=="C" && Cprev==-1{
				Cprev=num
				}
			if at.Name=="N" && Cprev!=-1 && N==-1 && at.Molid>M.Atoms[Cprev].Molid{
				N=num
				}
			if at.Name=="C" && Cprev!=-1 && at.Molid>M.Atoms[Cprev].Molid{
				C=num
				}
			if at.Name=="CA" && Cprev!=-1 && at.Molid>M.Atoms[Cprev].Molid{
				Ca=num
				}
			if at.Name=="N" && Ca!=-1 && at.Molid>M.Atoms[Ca].Molid{
				Npost=num
				}
			//when we have them all, save an unit
			if Cprev!=-1 && Ca!=-1 && N!=-1 && C !=-1 && Npost!=-1 {
				//We check that the residue ids are what they are supposed to be
				r1:=M.Atoms[Cprev].Molid
				r2:=M.Atoms[N].Molid
				r2a:=M.Atoms[Ca].Molid
				r2b:=M.Atoms[C].Molid
				r3:=M.Atoms[Npost].Molid
				if r1!=r2-1 || r2!=r2a || r2a!=r2b || r2b != r3-1{
					return  nil, fmt.Errorf("Incorrect backbone")
					}
				if (len(resran)==2 && (r2>=resran[0] && r2<=resran[1])) || IsInInt(resran,r2){
					temp:=[]int{Cprev,N,Ca,C,Npost}
					RamaList=append(RamaList,temp)
					}
				N=Npost
				Ca=-1
				Cprev=C
				C=-1
				Npost=-1					
				}
			}
		}
//	fmt.Println("Rama",Rama, "failed", failed)
	return RamaList,nil	
	}
