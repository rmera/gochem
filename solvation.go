/*
 * solvation.go, part of gochem
 *
 * Copyright 2020 Raul Mera A. (raulpuntomeraatusachpuntocl)
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
 *
*/
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/
package chem

import (
	"fmt"
	"log"
	"math"
	"runtime"
	"sort"
	"strings"

	//	"sort"
	//	"strconv"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

//Calculates the RDF for a trajectory given the indexes of the solute atoms, the solvent molecule name, the step for the "layers" and the cutoff.
//The code includes some extra comments, so it can be used as a template for concurrent trajectory processing.
func ConcMolRDF(traj ConcTraj, mol *Molecule, refindexes []int, residues []string, step, end float64, com ...bool) ([]float64, []float64, error) {
	var ret []float64
	simulframes := runtime.NumCPU()
	runtime.GOMAXPROCS(simulframes) //shouldn't be needed anymore, I think.
	frames := make([]*v3.Matrix, simulframes, simulframes)
	framesread := 0
	for i, _ := range frames {
		frames[i] = v3.Zeros(traj.Len())
	}
	results := make([]chan []float64, len(frames))
	for i, _ := range results {
		results[i] = make(chan []float64)
	}
	var err error
	for i := 0; ; i++ {
		if err != nil { //if we got a LastFrameError in the previous
			break
		}
		coordchans, err := traj.NextConc(frames) //we get a slice of channels, each of which will receive a frame. They are sorted by the frame they receive.
		if err != nil {
			if err, ok := err.(LastFrameError); ok {
				if coordchans == nil {
					break
				}
			} else {
				if err, ok := err.(Error); ok {
					err.Decorate(fmt.Sprintf("ConcMolRDF, failed when reading the %d th frame", i))
					return nil, nil, err
				}
				return nil, nil, err // somehow it wasn't a chem.Error.  This should never happen.
			}
		}
		for key, channel := range coordchans {
			go unitRDF(channel, results[key], mol, refindexes, residues, step, end, com...) //we give each of the channels we got from NextConc to a gorutine "worker" that performs the analysis. We pass them a chan to transmit the results back.
		}
		//Here we go through the "results" channels, which are sorted by frame, so if we just iterate with a for, we'll get the results for the frames in the right order (not that it matters here).
		for _, k := range results {
			if k == nil {
				break //shouldn't happen
			}
			rettemp := <-k
			if len(rettemp) == 0 { //we ran out of frames
				break
			}
			if ret == nil {
				ret = make([]float64, len(rettemp))
			}
			for i, v := range ret {
				ret[i] = v + rettemp[i]
			}
			framesread++
		}
	}
	ret, ret2 := mdfFromcdf(ret, framesread, step)
	return ret, ret2, nil
}

//The worker function for the RDF
func unitRDF(channelin chan *v3.Matrix, channelout chan []float64, mol Atomer, refindexes []int, residues []string, step, end float64, com ...bool) {
	if channelin != nil {
		temp := <-channelin
		rdf := FrameUMolCRDF(temp, mol, refindexes, residues, step, end, com...)
		channelout <- rdf
	} else {
		channelout <- nil
	}
	return
}

func MolRDF(traj Traj, mol Atomer, refindexes []int, residues []string, step, end float64, frameskip int, com ...bool) ([]float64, []float64, error) {
	var ret []float64
	coords := v3.Zeros(mol.Len())
	framesread := 0
	var err error
	if frameskip < 1 {
		frameskip = 1
	}
reading:
	for i := 0; ; i++ {
		if i > 0 && i%frameskip != 0 && err == nil {
			err = traj.Next(nil) //if this err is not nil, the next traj.Next() will not be excecuted, whether it's a skip or a read. Instead, we'll go directly to error processing.
			continue
		} else if err == nil { //in case the frame we skipped before gave an error
			err = traj.Next(coords)
		}
		if err != nil {
			switch err := err.(type) {
			case LastFrameError:
				break reading
			case Error:
				err.Decorate(fmt.Sprintf("MolRDF: Failed while reading the %d th frame", i))
				return nil, nil, err
			default:
				return nil, nil, err

			}
		}
		rdf := FrameUMolCRDF(coords, mol, refindexes, residues, step, end, com...) ///only difference
		if ret == nil {
			ret = make([]float64, len(rdf))
		}
		for j, _ := range ret {
			ret[j] = ret[j] + rdf[j]
		}
		framesread++

	}
	ret, ret2 := mdfFromcdf(ret, framesread, step)

	return ret, ret2, nil
}

// transforms a cdf into mdf, and normalizes it by the total value
//it also divides from the total number of datapoints to go from
//accumulated data, to the mean. It returns a slice with the molecule
//density (divide by the volume of the sector) and a slice with the
//number of molecules i nthe sector.
// it overwrites the original slice!
func mdfFromcdf(ret []float64, framesread int, step float64) ([]float64, []float64) {
	ret2 := make([]float64, len(ret))
	vp := (4.0 / 3.0) * math.Pi
	//vp := 4 * math.Pi
	//	avtotalsolv := ret[len(ret)-1]   //I think this is not needed, as we are already dealing with densities
	var vol float64
	acc := 0.0
	for i, _ := range ret {
		if i >= 1 {
			acc = acc + ret[i-1]
			ret2[i] = ret[i]
			ret[i] = ret[i] - acc
		}
	}
	for i, _ := range ret {
		if i >= 1 {
			fi := float64(i)
			vol = vp * (math.Pow((fi+1)*step, 3) - math.Pow((fi)*step, 3))

		} else if i == 0 {
			vol = vp * math.Pow(step, 3)

		} else {
			vol = 1
		}
		ret[i] = ret[i] / float64(framesread)
		ret2[i] = ret2[i] / float64(framesread)
		ret[i] = ret[i] / vol

	}
	last := ret[len(ret)-1]
	for i, _ := range ret {
		ret[i] = ret[i] / last
	}

	return ret, ret2
}

//Obtains the Unnormalized "Cummulative Molecular RDF" for one solvated structure. The RDF would be these values averaged over several structures.
func FrameUMolCRDF(coord *v3.Matrix, mol Atomer, refindexes []int, residues []string, step, end float64, com ...bool) []float64 {
	if step <= 0 {
		step = 0.1
	}
	if end <= 0 {
		end = 15
	}
	totalsteps := int(end / step)
	ret := make([]float64, 0, totalsteps)
	// Here we get the RDF by counting the molecules at X distance or less from the solute, and subtracting the result from the previous callculation.
	//This means that we do more calculation than needed, as every time keep including the solvent that was in previous layers.
	//	acclen := 0.0

	res := DistRank(coord, mol, refindexes, residues, end, com...)
	sort.Sort(res)
	dists := res.Distances()
	for i := 1; i <= totalsteps; i++ {
		limit := float64(i) * step
		n := 0
		for _, v := range dists {
			if v > limit {
				break
			}
			n++

		}
		ret = append(ret, float64(n))

	}
	return ret
}

type resAndChain struct {
	ResID int
	Chain string
}

//returns all the residue numbers in mol covered by indexes
func allResIDandChains(mol Atomer, indexes []int) []*resAndChain {
	ret := make([]*resAndChain, 0, 2)
	for _, i := range indexes {
		at := mol.Atom(i)
		if !repeated(at.MolID, at.Chain, ret) {
			ret = append(ret, &resAndChain{ResID: at.MolID, Chain: at.Chain})
		}
	}
	return ret
}

func repeated(id int, chain string, rac []*resAndChain) bool {
	for _, v := range rac {
		if v.Chain == chain && id == v.ResID {
			return true
		}
	}
	return false
}

//DistRank determines, for a reference set of coordinates and a set of residue names, the minimum distance between any atom from the reference
//and any atom from each residue with one of the names given (or the centroid of each residue, if a variadic "com" bool is given)
//returns a list with ID and distances, which satisfies the sort interface and has several other useful methods.
func DistRank(coord *v3.Matrix, mol Atomer, refindexes []int, residues []string, cutoff float64, com ...bool) MolDistList {
	ranks := make([]*molDist, 0, 30)
	//	resIDs := make([]*resAndChain, 0, 30)
	var molname string
	var id, molid_skip int //The "skip" variables keep the residue just being read or discarded to avoid reading a residue twice.
	molid_skip = -1
	var chain, chain_skip string
	var distance float64
	water := v3.Zeros(3)
	var ref *v3.Matrix
	ownresIDs := allResIDandChains(mol, refindexes) //We could call this upstream and just get this numbers, but I suspect
	ref = v3.Zeros(len(refindexes))
	ref.SomeVecs(coord, refindexes)
	if cutoff < 0 {
		cutoff = 10 //if a negative cutoff is given we use 10 A as the default.
	}
	cutoffplus := cutoff + 3
	//	chunk := NewTopology(0, 1)
	tmp := v3.Zeros(1)
	//	fmt.Println("Start looking!") ///////////////////
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		molname = at.Molname
		id = at.MolID
		chain = at.Chain
		var test *v3.Matrix
		//a little pre-screening for waters
		if isInString([]string{"SOL", "WAT", "HOH"}, molname) && dist(ref.VecView(0), coord.VecView(i), tmp) > cutoffplus {
			continue

		}
		if isInString(residues, molname) && !repeated(id, chain, ownresIDs) && (id != molid_skip || chain != chain_skip) {
			expectedreslen := 6
			indexes := make([]int, 1, expectedreslen)
			indexes[0] = i
			//	fmt.Println("New mol!", molname, id, chain) //////////
			for j := i + 1; j < mol.Len(); j++ {
				at2 := mol.Atom(j)
				if at2.MolID != id {
					break
				}
				indexes = append(indexes, j)
			}
			if len(indexes) == 3 {
				test = water
			} else {
				test = v3.Zeros(len(indexes))
			}
			test.SomeVecs(coord, indexes)
			//This is a bit ugly
			if len(com) != 0 && com[0] == true {
				var mass *mat.Dense
				topol := NewTopology(0, 1)
				topol.SomeAtoms(mol, indexes)
				massp, err := topol.Masses()
				if err != nil {
					mass = mat.NewDense(len(massp), 1, massp)
				}
				c, err := CenterOfMass(test, mass)
				if err != nil { //this really is very unlikely to fail. The ref matrix would have to be wrong. It could even deserve a panic.
					com[0] = false //if it fails, we don't try again, for consistency. Any previous successful COM use will not comparable with numbers used from now on.
					log.Println("Couldn't obtain the COM/centroid. Worked with all solvent atoms")
					distance = MolShortestDist(test, ref) //we don't panic, we just keep going with all atoms
				}
				distance = MolShortestDist(c, ref)
			} else {
				distance = MolShortestDist(test, ref)
			}
			if distance <= cutoff {
				ranks = append(ranks, &molDist{Distance: distance, MolID: id})
			}
			molid_skip = id
			chain_skip = chain

		}
	}
	return MolDistList(ranks)

}

type molDist struct {
	Distance float64
	MolID    int
}

func (M *molDist) str() string {
	return fmt.Sprintf("D: %4.3f ID: %d", M.Distance, M.MolID)
}

type MolDistList []*molDist

func (M MolDistList) Swap(i, j int) {
	M[i], M[j] = M[j], M[i]
}
func (M MolDistList) Less(i, j int) bool {
	return M[i].Distance < M[j].Distance
}
func (M MolDistList) Len() int {
	return len(M)
}

func (M MolDistList) Distance(i int) float64 {
	return M[i].Distance
}

func (M MolDistList) MolID(i int) int {
	return M[i].MolID
}

func (M MolDistList) String() string {
	retslice := make([]string, len(M))
	for i := range M {
		retslice[i] = M[i].str()
	}
	return strings.Join(retslice, "\n")
}

//Returns a list of the molIDs in al the lists
func (M MolDistList) Distances() []float64 {
	ret := make([]float64, len(M))
	for i := range M {
		ret[i] = M[i].Distance //This absolutely should never fail!
	}
	return ret
}

//Returns a list of the molIDs in al the lists
func (M MolDistList) MolIDs() []int {
	ret := make([]int, len(M))
	for i := range M {
		ret[i] = M[i].MolID //This absolutely should never fail!
	}
	return ret
}

//Returns a list with all the MolIDs and a list with all the distances
func (M MolDistList) Data() ([]int, []float64) {
	ret := make([]int, len(M))
	ret2 := make([]float64, len(M))
	for i := range M {
		ret[i] = M[i].MolID //This absolutely should never fail!
		ret2[i] = M[i].Distance
	}
	return ret, ret2
}

//Returns a list of all the atom IDs for all the residues in the list.
//Only a convenience function.
func (M MolDistList) AtomIDs(mol Atomer) []int {
	molids := M.MolIDs()
	return Molecules2Atoms(mol, molids, []string{})
}

//Aggregates the receiver and the arguments in the receiver
//it removes entries with repeated MolIDs
func (M MolDistList) Merge(list ...MolDistList) {
	for _, v := range list {
		ref := M.MolIDs()
		for _, w := range v {
			if !isInInt(ref, w.MolID) {
				M = append(M, w)
			}
		}
	}

}

//MolShortestDistGiven two sets of coordinates, it obtains the shortest distance from any 2 points
//in the set. This is probably not a very efficient way to do it.
func MolShortestDist(test, ref *v3.Matrix) float64 {
	temp := v3.Zeros(1)
	var d1, dclosest float64
	var vt1, vr1 *v3.Matrix // vtclosest,vr1, vrclosest *v3.Matrix
	dclosest = 100000       //we initialize this to some crazy large value so it's immediately replaced with the first calculated distance.
	// dist(ref.VecView(0), test.VecView(0), temp)
	for i := 0; i < test.NVecs(); i++ {
		vt1 = test.VecView(i)
		for j := 0; j < ref.NVecs(); j++ {
			vr1 = ref.VecView(j)
			temp.Sub(vr1, vt1)
			d1 = temp.Norm(2)
			if d1 < dclosest {
				dclosest = d1
			}
		}
	}
	return dclosest
}

func dist(r, t, temp *v3.Matrix) float64 {
	temp.Sub(r, t)
	return temp.Norm(2)
}
