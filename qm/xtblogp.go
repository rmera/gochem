/*
 * xtblogp.go, part of gochem.
 *
 *
 * Copyright 2016 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 *
 */
//In order to use this part of the library you need the xtb program, which must be obtained from Prof. Stefan Grimme's group.
//Please cite the the xtb references if you used the program.

package qm

import (
	"fmt"
	"log"
	"math"
	"os"
	"slices"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

const (
	defgmethod   = "gfn2"
	defmdmethod  = "gfnff"
	defmdtime    = 5000
	defframeskip = 1000
	defimgthres  = 0.1
	defmdtemp    = 310.0
	defminn      = 5

//	defmaxrelsigma = 0.1 //max relative stdev for each deltaG average
)

var methodsav = []string{"gfn0", "gfn2", "gfnff"}

type LPOptions struct {
	MinN      int    //Minimum N frames for each G (water and octanol). The value for logP will still be returned but with an error.
	GMethod   string //gfn2, gfn0 and gfnff
	MDMethod  string //same
	MDTime    int    //0 means no MD
	MDTemp    float64
	FrameSkip int //every how many steps to sample. 0 is 1 sample at the end

	ImaginaryFreqThres float64 //how large can an imaginary mode be before we decide the structure is not a minimum?
}

func (O *LPOptions) Check() {
	if !slices.Contains(methodsav, O.GMethod) {
		log.Printf("Selected G method %s not available, will use the default %s\n", O.GMethod, defgmethod)
		O.GMethod = defgmethod
	}

	if !slices.Contains(methodsav, O.MDMethod) {
		log.Printf("Selected MD method %s not available, will use the default %s\n", O.MDMethod, defmdmethod)
		O.MDMethod = defmdmethod
	}
	if O.MDTime <= 0 {
		log.Printf("Invalid MD time selected %d. Will use the default: %d", O.MDTime, defmdtime)
		O.MDTime = defmdtime
	}
	if O.FrameSkip <= 0 {
		log.Printf("Invalid skip frame rate selected %d. Will use the default: %d", O.FrameSkip, defframeskip)
		O.FrameSkip = defframeskip
	}
	if O.MDTemp <= 0 {
		log.Printf("Invalid MD temperature selected %5.1f. Will use the default: %5.1f", O.MDTemp, defmdtemp)
		O.FrameSkip = defframeskip
	}

	if O.MDTemp <= 0 {
		log.Printf("Invalid MD temperature selected %5.1f. Will use the default: %5.1f", O.MDTemp, defmdtemp)
		O.FrameSkip = defframeskip
	}
	if O.MinN <= 0 {
		log.Printf("Invalid Minimum N for each deltaG average %d. Will use the default: %d", O.MinN, defminn)
		O.FrameSkip = defframeskip
	}

	if O.ImaginaryFreqThres <= 0 {
		log.Printf("Invalid imaginary frequency threshold for minumum %5.2f. Will use the default: %5.2f", O.ImaginaryFreqThres, defimgthres)
		O.ImaginaryFreqThres = defimgthres
	}
}

func (O *LPOptions) SetDefaults() {
	O.GMethod = defgmethod
	O.MDMethod = defmdmethod
	O.MDTime = defmdtime
	O.FrameSkip = defframeskip
	O.MDTemp = defmdtemp
	O.MinN = defminn
	O.ImaginaryFreqThres = defimgthres
}

func LogP(coords *v3.Matrix, mol chem.AtomMultiCharger, options ...*LPOptions) (float64, error) {
	var o *LPOptions
	if len(options) == 0 {
		o = new(LPOptions)
		o.SetDefaults()
	} else {
		o = options[0]
		o.Check()
	}
	if mol.Charge() != 0 {
		log.Printf("LogP is defined for non-ionized molecules!. Will procede regardless\n")
	}
	Q := new(Calc)
	Q.Job = &Job{MD: true}
	Q.Method = o.MDMethod
	Q.MDTime = o.MDTime
	Q.MDTemp = o.MDTemp
	watere := 80.0
	octanole := 1000.0 //this is a rather poor hack in goChem to have wet octanol and octanol at the same time. Wet octanol
	//just got assigned the "dielectric" 1000.
	solvents := []float64{watere, octanole}
	Gs := []float64{}
	var lackofsamples error
	for _, v := range solvents {
		Q.Dielectric = v
		solvname := "water"
		if v > 900.0 {
			solvname = "n-octanol(wet)"
		}
		xtb := NewXTBHandle()
		err := xtb.BuildInput(coords, mol, Q)
		if err != nil {
			return 0, fmt.Errorf("couldn't build xtb input: %s", err.Error())
		}
		err = xtb.Run(true)
		if err != nil {
			return 0, fmt.Errorf("xtb didn't run correctly: %s", err.Error())
		}
		Q.Method = o.GMethod //sorry, it's late
		gs, err := xtbAvG(o.FrameSkip, Q, o.ImaginaryFreqThres)
		Q.Method = o.MDMethod
		if err != nil {
			return 0, fmt.Errorf("Couldn't obtain free energy for solvent %s: %v", solvname, err)
		}
		if len(gs) < o.MinN {
			es := fmt.Sprintf("Not enough samples for deltaG in %s. Wanted: %d. Got: %d\n", solvname, o.MinN, len(gs))
			if lackofsamples == nil {
				lackofsamples = fmt.Errorf(es)
			} else {
				lackofsamples = fmt.Errorf("%v %s", lackofsamples, es)
			}
		}
		Gs = append(Gs, boltzmannAv(gs, o.MDTemp))
	}
	if len(Gs) != 2 {
		//This has to be a bug. I'll panic so I get a bug report if it happens.
		fs := fmt.Sprintf
		panic(fs("goChem/LogP: For some reason I don't have both, or have too many (%d) of the deltaGs I need for the LogP calculation", len(Gs)))
	}
	deltaG := Gs[1] - Gs[0]
	P := math.Exp(-deltaG / (chem.R * o.MDTemp))
	return math.Log10(P), lackofsamples

}

func xtbAvG(skip int, Q *Calc, imgthres ...float64) ([]float64, error) {
	ithres := 0.1
	if len(imgthres) > 0 {
		ithres = imgthres[0]
	}
	mol, trj, err := chem.XYZFileAsTraj("xtb.trj")
	if err != nil {
		return nil, fmt.Errorf("Couldn't open xtb trajectory: %s", err.Error())
	}
	cwd, err := os.Getwd()
	if err != nil {
		panic(fmt.Sprintf("goChem/LogP/xtbAvG: Couldn't identify the current directory!: %s", err.Error()))
	}
	//The following 2 are to ensure that we return the environment to the way it was when we got it.
	defer os.Chdir(cwd)
	defer func(Q *Calc) { Q.Job = &Job{MD: true} }(Q)

	Q.Job = &Job{Forces: true}
	coord := v3.Zeros(trj.Len())
	Gs := make([]float64, 5)
	totGsq := 0.0
	for i := 0; ; i++ {
		err = trj.Next(coord)
		if err != nil {
			if err.Error() == "EOF" {
				break
			}
			return nil, fmt.Errorf("Failed to read (and, possibly, discard), frame %d: %s", i, err.Error())
		}
		if i%skip != 0 {
			continue
		}

		dirname := fmt.Sprintf("%d", i)
		xtb := NewXTBHandle()
		os.Mkdir(dirname, 0777)
		os.Chdir(dirname)
		xtb.SetName(fmt.Sprintf("G_%d", i))
		xtb.BuildInput(coord, mol, Q)
		err = xtb.Run(true)
		if err != nil {
			log.Printf("Failed to run Forces calculation for frame %d: %v\n", i, err)
			continue
		}
		_, err = xtb.OptimizedGeometry(mol)
		if err != nil {
			log.Printf("Failed to obtain optimized geometry for frame %d: %v", i, err)
			continue
		}
		img, err := xtb.LargestImaginary()
		if err != nil {
			log.Printf("Skipping frame %d because of error in reading frequencies: %v", i, err)
			continue
		}
		if img > ithres {
			log.Printf("Skipping frame %d because of negative frequency. Wavenumber: %5.3f", i, img*-1)
			os.Rename("vibspectrum", fmt.Sprintf("vibspectrum_%d", i))
			continue
		}
		G, err := xtb.FreeEnergy()
		if err != nil {
			log.Printf("Failed to obtain free energy for frame %d: %v", i, err)
			continue
		}
		fmt.Printf("G: %5.3f Frame: %d\n", G, i)
		Gs = append(Gs, G)
		totGsq += (G * G)
		os.Chdir(cwd)
	}

	return Gs, nil

}

func boltzmannAv(es []float64, T float64) float64 {
	p := func(e float64) float64 { return math.Exp(-e / (chem.R * T)) }
	ps := make([]float64, 0, len(es))
	Z := 0.0
	for _, v := range es {
		pi := p(v)
		ps = append(ps, pi)
		Z += pi
	}
	bav := 0.0
	for i, v := range ps {
		bav += es[i] * (v / Z)
	}
	return bav

}
