/*
 * plot.go, part of gochem.
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



import (
		"fmt"
//		"github.com/skelterjohn/go.matrix"
  //      "code.google.com/p/plotinum/vg"
        "code.google.com/p/plotinum/plot"
        "code.google.com/p/plotinum/plotter"
        "image/color"
)

func RamachandranPlot(fulldata [][][]float64) error{
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
		p.Title.Text = "Ramachandran plot"
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
		filename:=fmt.Sprintf("test/Ramachandran%d.png",number)
		if err := p.Save(4, 4, filename); err != nil {
			return err
		}
	}
	return nil
}

