// +build plot

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
	"code.google.com/p/plotinum/plot"
	"code.google.com/p/plotinum/plotter"
	"fmt"
	"image/color"
	"math"
	"strings"
)

type RamaSet struct {
	Cprev int
	N int
	Ca int
	C int
	Npost int
	Molid int
	Molname string
	}

/*Produce plots, in png format for the ramachandran data (psi and phi dihedrals)
  contained in fulldata, which can contain data for various different snapshopts.
  In the latter case, many png files are produced. The file names are plotnameXX.png
  where XX is the frame number (not limited to digits). Returns an error*/
func RamaPlot(data [][]float64, plotname, title string, tag int) error  {
	if data == nil {
		panic("Given nil data")
	}
/*	pts := make(plotter.XYs, len(data)) //just first frame for now
	//this might not be too efficient
	for key, val := range data {
		pts[key].X = val[0]
		pts[key].Y = val[1]
	}
*/
	// Create a new plot, set its title and
	// axis labels.
	p, err := plot.New()
	if err != nil {
		return err
	}
	p.Title.Text = title //"Ramachandran plot"
	p.X.Label.Text = "Phi"
	p.Y.Label.Text = "Psi"
	//Constant axes
	p.X.Min = -180
	p.X.Max = 180
	p.Y.Min = -180
	p.Y.Max = 180
	// Draw the grid
	p.Add(plotter.NewGrid())
	critical:=0 // this is the residue where I run out of RB combinations and have to start adding green.
	temp:=make(plotter.XYs,1)
	//Here we try to produce one color for each point. First e go from red to blue, and we continue with blue to green.
	for key,val:=range(data){
		temp[0].X=val[0]
		temp[0].Y=val[1]
		// Make a scatter plotter and set its style.
		s, err := plotter.NewScatter(temp) //(pts)
		if err != nil {
			return err
		}
		var g uint8
		norm:=(2*255.0/len(data))+1
		b:=uint8(key*norm)
		r:=uint8(255)-b
		if norm==256{
			critical=key
			}
		if key*norm>255{
			g=uint8(norm*(key-critical))
			b=255-g
			r=0
		}
		if key==tag{
			s.GlyphStyle.Shape=plot.PyramidGlyph{}
			}
		s.GlyphStyle.Color = color.RGBA{R:r,B:b,G:g,A: 255}
//		fmt.Println(r,b,g, key, norm, len(data)) //////////////////////////
		// Add the plotter
		p.Add(s)
	}
	// Save the plot to a PNG file.
	filename := fmt.Sprintf("%s.png", plotname)
	if err := p.Save(4, 4, filename); err != nil {
		return err
	}
	return nil
}

/*Obtain psi and phi angles for the molecule M, considerint the dihedrals in [][]int
  for all the frames in frames. It returns a slice of slices (one per frame) of slices
  (with 2 elemtns), for psi and phi angles, and an error*/
func RamaCalc(M *CoordMatrix, dihedrals []RamaSet) ([][]float64, error) {
	if M == nil || dihedrals == nil {
		return nil, fmt.Errorf("RamaCalc: Given nil data")
	}
	r,_:=M.Dims()
	Rama := make([][]float64, 0, len(dihedrals))
	for _, j := range dihedrals {
		if j.Npost>=r{
			return nil, fmt.Errorf("RamaCalc: Index requested out of range")
		}
		Cprev := RowView(M, j.Cprev)
		N := RowView(M, j.N)
		Ca := RowView(M, j.Ca)
		C := RowView(M, j.C)
		Npost := RowView(M, j.Npost)
		phi := Dihedral(Cprev, N, Ca, C)
		psi := Dihedral(N, Ca, C, Npost)
		temp := []float64{phi * (180 / math.Pi), psi * (180 / math.Pi)}
		Rama = append(Rama, temp)
	}
	return Rama, nil
}

//Filter the set of dihedral angles of a ramachandran plot by residue.(ex. only GLY, everything but GLY)
//The 3 letter code of the residues to be filtered in or out is in filterdata, whether they are filter in 
//or out depends on shouldBePresent. It returns the filtered data and a slice containing the indexes in 
//the new data of the residues in the old data, when they are included, or -1 when they are not included.
func RamaResidueFilter(dihedrals []RamaSet, filterdata []string, shouldBePresent bool) ([]RamaSet, []int){
	RetList := make([]RamaSet, 0, 0)
	Index:=make([]int,len(dihedrals))
	var added int
	for key,val:=range(dihedrals){
		isPresent:=isInString(filterdata, val.Molname)
		if isPresent==shouldBePresent{
			RetList=append(RetList,val)
			Index[key]=added
			added++
			}else{
			Index[key]=-1
			}
		}
	return RetList, Index
}


/*RamaList takes a molecule and obtains a list of lists of five int. Each element
  contain the indexes needed for one dihedral of a Rama plot. It gets the dihedral
  indices for all residues in the range resran, if resran has 2 elements defining the
  boundaries. Otherwise, returns dihedral lists for the residues included in
  resran. If resran has 2 elements and the last is -1, RamaList will
  get all the dihedral for residues from resran[0] to the end of the chain.
  It only obtain dihedral lists for residues belonging to a chain included in chains */
func RamaList(M Ref, chains string, resran []int) ([]RamaSet, error) {
	RamaList := make([]RamaSet, 0, 0)
	if len(resran) == 2 {
		if resran[1] == -1 {
			resran[1] = 999999999 //should work!
		}
	}
	if M == nil {
		return nil, fmt.Errorf("nil Molecule")
	}
	C := -1
	N := -1
	Ca := -1
	Cprev := -1
	Npost := -1
	chainprev := byte('9') //any non-valid chain name
	for num := 0; num < M.Len(); num++ {
		at := M.Atom(num)
		//First get the indexes we need
		if strings.Contains(chains, string(at.Chain)) || at.Chain == ' ' {
			if at.Chain != chainprev {
				chainprev = at.Chain
				C = -1
				N = -1
				Ca = -1
				Cprev = -1
				Npost = -1
			}
			if at.Name == "C" && Cprev == -1 {
				Cprev = num
			}
			if at.Name == "N" && Cprev != -1 && N == -1 && at.Molid > M.Atom(Cprev).Molid {
				N = num
			}
			if at.Name == "C" && Cprev != -1 && at.Molid > M.Atom(Cprev).Molid {
				C = num
			}
			if at.Name == "CA" && Cprev != -1 && at.Molid > M.Atom(Cprev).Molid {
				Ca = num
			}
			if at.Name == "N" && Ca != -1 && at.Molid > M.Atom(Ca).Molid {
				Npost = num
			}
			//when we have them all, we save
			if Cprev != -1 && Ca != -1 && N != -1 && C != -1 && Npost != -1 {
				//We check that the residue ids are what they are supposed to be
				r1 := M.Atom(Cprev).Molid
				r2 := M.Atom(N).Molid
				r2a := M.Atom(Ca).Molid
				r2b := M.Atom(C).Molid
				r3 := M.Atom(Npost).Molid
				if r1 != r2-1 || r2 != r2a || r2a != r2b || r2b != r3-1 {
					return nil, fmt.Errorf("Incorrect backbone")
				}
				if (len(resran) == 2 && (r2 >= resran[0] && r2 <= resran[1])) || isInInt(resran, r2) {
					temp := RamaSet{Cprev, N, Ca, C, Npost, r2, M.Atom(Ca).Molname}
					RamaList = append(RamaList, temp)
				}
				N = Npost
				Ca = -1
				Cprev = C
				C = -1
				Npost = -1
			}
		}
	}
	//	fmt.Println("Rama",Rama, "failed", failed)
	return RamaList, nil
}
