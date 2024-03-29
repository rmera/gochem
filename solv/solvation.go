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

package solv

import (
	"fmt"
	"log"
	"math"
	"runtime"
	"sort"
	"strings"

	//	"sort"
	//	"strconv"
	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

const epsilon = 0.0001

func EllipsoidAxes(coords *v3.Matrix, epsilon float64, mol ...chem.Masser) ([]float64, error) {
	var masses []float64
	var err2 error
	var err error
	if len(mol) < 0 {
		masses = nil
	} else {
		masses, err = mol[0].Masses()
		if err != nil {
			masses = nil
			err2 = err
		}
	}
	moment, err := chem.MomentTensor(coords, masses)
	if err != nil {
		return nil, err
	}
	rhos, err := chem.Rhos(moment, epsilon)
	if err != nil {
		return nil, err
	}

	return rhos, err2
}

// Options contains options for the RDF/MDDF calculation
type Options struct {
	com  bool
	cpus int
	step float64
	end  float64
	skip int
}

// Returns a Options with the default options.
func DefaultOptions() *Options {
	ret := new(Options)
	ret.com = false
	ret.cpus = runtime.NumCPU()
	ret.step = 0.1
	ret.end = 10
	ret.skip = 1
	return ret
}

// Returns whether to use center of mass for solvent in the calculations
// and sets the value to the one given, if any
func (r *Options) COM(com ...bool) bool {
	ret := r.com
	if len(com) > 0 {
		r.com = com[0]
	}
	return ret
}

// Returns the current value of the Cpus options (the number of gorutines to
// use on the concurrent calculation) and sets it, if
// a valid value is given
func (r *Options) Cpus(cpus ...int) int {
	ret := r.cpus
	if len(cpus) > 0 && cpus[0] > 0 {
		r.cpus = cpus[0]
	}
	return ret
}

// Returns the skipped frames (for functions where it's applicable) and sets it, if
// a valid value is given
func (r *Options) Skip(skip ...int) int {
	ret := r.skip
	if len(skip) > 0 && skip[0] > 0 {
		r.skip = skip[0]
	}
	return ret
}

// Returns the distance step to be used in the RDF/MDDF calculation
// and sets if to a value, if a valid value is given
func (r *Options) Step(step ...float64) float64 {
	ret := r.step
	if len(step) > 0 && step[0] > 0 {
		r.step = step[0]
	}
	return ret
}

// Returns the maximum distance from the solute to be considered
// in the RDF/MDDF calculation and sets if to a value, if given
func (r *Options) End(end ...float64) float64 {
	ret := r.end
	if len(end) > 0 && end[0] > 0 {
		r.end = end[0]
	}
	return ret
}

// ConcMolRDF calculates the RDF for a trajectory given the indexes of the solute atoms, the solvent molecule name, the step for the "layers" and the cutoff.
// It processes several frames of the trajectory concurrently, depending on the logical CPUs available.
// The code includes some extra comments, so it can be used as a template for concurrent trajectory processing.
func ConcMolRDF(traj chem.ConcTraj, mol *chem.Molecule, refindexes []int, residues []string, options ...*Options) ([]float64, []float64, error) {
	var o *Options
	if len(options) > 0 {
		o = options[0]
	} else {
		o = DefaultOptions()
	}
	A := 1.0
	B := 1.0
	if mol.Len() > 1 {
		scoords := v3.Zeros(len(refindexes))
		scoords.SomeVecs(mol.Coords[0], refindexes)
		smol := chem.NewTopology(0, 1)
		smol.SomeAtoms(mol, refindexes)
		elip, err := EllipsoidAxes(scoords, epsilon, smol)
		if err != nil {
			return nil, nil, err
		}
		//fmt.Println(elip, elip[0]*elip[1]*elip[2]) /////////////////////////////////
		A = elip[0] / elip[2]
		B = elip[1] / elip[2]

	}

	var ret []float64
	frames := make([]*v3.Matrix, o.cpus, o.cpus)
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
			if err, ok := err.(chem.LastFrameError); ok {
				if coordchans == nil {
					break
				}
			} else {
				if err, ok := err.(chem.Error); ok {
					err.Decorate(fmt.Sprintf("ConcMolRDF, failed when reading the %d th frame", i))
					return nil, nil, err
				}
				return nil, nil, err // somehow it wasn't a chem.Error.  This should never happen.
			}
		}
		for key, channel := range coordchans {
			go unitRDF(channel, results[key], mol, refindexes, residues, o) //we give each of the channels we got from NextConc to a gorutine "worker" that performs the analysis. We pass them a chan to transmit the results back.
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
	//fmt.Println(A, B) /////////////
	ret, ret2, err := MDFFromCDF(ret, framesread, A, B, o.step)
	return ret, ret2, err
}

// The worker function for the RDF
func unitRDF(channelin chan *v3.Matrix, channelout chan []float64, mol chem.Atomer, refindexes []int, residues []string, o *Options) {
	if channelin != nil {
		temp := <-channelin
		rdf := FrameUMolCRDF(temp, mol, refindexes, residues, o)
		channelout <- rdf
	} else {
		channelout <- nil
	}
	return
}

// MolRDF calculates the RDF for a trajectory given the indexes of the solute atoms, the solvent molecule name, the step for the "layers" and the cutoff.
// API BREAK: mol used to be chem.Atomer.
func MolRDF(traj chem.Traj, mol *chem.Molecule, refindexes []int, residues []string, options ...*Options) ([]float64, []float64, error) {
	var o *Options
	if len(options) > 0 {
		o = options[0]
	} else {
		o = DefaultOptions()
	}
	var ret []float64
	coords := v3.Zeros(mol.Len())
	framesread := 0
	var err error
	A := 1.0
	B := 1.0

reading:
	for i := 0; ; i++ {
		if i > 0 && i%o.skip != 0 && err == nil {
			err = traj.Next(nil) //if this err is not nil, the next traj.Next() will not be excecuted, whether it's a skip or a read. Instead, we'll go directly to error processing.
			continue
		} else if err == nil { //in case the frame we skipped before gave an error
			err = traj.Next(coords)
		}
		if err != nil {
			switch err := err.(type) {
			case chem.LastFrameError:
				break reading
			case chem.Error:
				err.Decorate(fmt.Sprintf("MolRDF: Failed while reading the %d th frame", i))
				return nil, nil, err
			default:
				return nil, nil, err

			}
		}
		rdf := FrameUMolCRDF(coords, mol, refindexes, residues, o) ///only difference
		if ret == nil {
			ret = make([]float64, len(rdf))
		}
		for j, _ := range ret {
			ret[j] = ret[j] + rdf[j]
		}
		if mol.Len() > 1 && framesread == 1 {
			elip, err := EllipsoidAxes(coords, epsilon, mol)
			if err != nil {
				return nil, nil, err
			}
			A = elip[0] / elip[2]
			B = elip[1] / elip[2]

		}

		framesread++

	}
	ret, ret2, err := MDFFromCDF(ret, framesread, A, B, o.step)

	return ret, ret2, err
}

// Obtains the radial standard-deviation distribution function from the unnormalized cummulative RDF and square RDF
func SQRDF2RSDF(avs, sqavs []float64, framesread int, step float64) []float64 {
	ret := make([]float64, len(avs))
	for i, v := range avs {
		ret[i] = sqavs[i]/float64(framesread) - math.Pow(v/float64(framesread), 2)

	}
	last := math.Sqrt(ret[len(ret)-1])
	for i, _ := range ret {
		ret[i] = math.Sqrt(ret[i]) / last
	}

	return ret

}

// MDDFFromCDF takes the sume of a cummulative distribution function over framesread frames.
// It obtaines the RDF/MDDF from there by dividing each "shell" by an approximation to the
// volume (the ellipsoid of inerta of the solute scaled to the corresponding radius) and divided
// by the frames read. It returns a slice with the average  density per shell (divided by volume)
// and another with the average number of molecules per shell, in both cases, divided by the value
// of the last shell. The function also requires the ratio of the largest semiaxes of the ellipsoid of inertia
// to its smallest semiaxis, A and B. It returns an error and nil slices if A and B are smaller than 1.
// it overwrites the original slice CDF slice!
// API BREAK: The original function did not take A and B, and simply approximated the volume as a sphere.
// This departs from the behavior described in the publication, in a way that could create wrong results.
// The current implementation also departs from the publication, which approximated the volume by a parallelepiped.
// I think this behavior is better, as it allows the MDDF to reduce to the RDF for spherically symmetric systems.
func MDFFromCDF(ret []float64, framesread int, A, B, step float64) ([]float64, []float64, error) {
	if A < 1.0 || B < 1.0 {
		return nil, nil, fmt.Errorf("goChem/solvation.MDFFromCDF: A and B are the ratio between the lartest axis of an elipsoid to the smallest, they can't be smaller than 1.0")
	}
	ret2 := make([]float64, len(ret))
	vp := (4.0 / 3.0) * math.Pi * A * B
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

	return ret, ret2, nil
}

// FrameUMolCRDF Obtains the the number of solvent molecules in each solvent shell, and the sqare of that number for each shell
// in a given frame.
func FrameUMolSQRDF(coord *v3.Matrix, mol chem.Atomer, refindexes []int, residues []string, options ...*Options) ([]float64, []float64) {
	//NOTE: It could be good to have this function take a slice of floats and put the results there, so as to avoid
	//allocating more than needed.
	var o *Options
	if len(options) > 0 {
		o = options[0]
	} else {
		o = DefaultOptions()
	}

	totalsteps := int(o.end / o.step)
	av := make([]float64, 0, totalsteps)
	sqav := make([]float64, 0, totalsteps)
	// Here we get the RDF by counting the molecules at X distance or less from the solute, and subtracting the result from the previous callculation.
	//This means that we do more calculation than needed, as every time keep including the solvent that was in previous layers.
	//	acclen := 0.0

	res := DistRank(coord, mol, refindexes, residues, o)
	sort.Sort(res)
	dists := res.Distances()
	var nprev float64 = 0
	for i := 1; i <= totalsteps; i++ {
		limit := float64(i) * o.step
		n := sort.SearchFloat64s(dists, limit)
		// the limit is not found, SearchFloat64 will return the index of the first element larger or equal than limit.
		// we want the one before that.
		if n > 0 {
			n--
		}
		av = append(av, float64(n)-nprev)
		sqav = append(av, math.Pow(float64(n)-nprev, 2))
		nprev = float64(n)

	}
	return av, sqav
}

// FrameUMolCRDF Obtains the Unnormalized Cummulative Molecular RDF for one solvated structure. The RDF would be these values averaged over several structures.
func FrameUMolCRDF(coord *v3.Matrix, mol chem.Atomer, refindexes []int, residues []string, options ...*Options) []float64 {
	//NOTE: It could be good to have this function take a slice of floats and put the results there, so as to avoid
	//allocating more than needed.
	var o *Options
	if len(options) > 0 {
		o = options[0]
	} else {
		o = DefaultOptions()
	}

	totalsteps := int(o.end / o.step)
	ret := make([]float64, 0, totalsteps)
	// Here we get the RDF by counting the molecules at X distance or less from the solute, and subtracting the result from the previous callculation.
	//This means that we do more calculation than needed, as every time keep including the solvent that was in previous layers.
	//	acclen := 0.0

	res := DistRank(coord, mol, refindexes, residues, o)
	sort.Sort(res)
	dists := res.Distances()
	for i := 1; i <= totalsteps; i++ {
		limit := float64(i) * o.step
		n := sort.SearchFloat64s(dists, limit)
		// the limit is not found, SearchFloat64 will return the index of the first element larger or equal than limit.
		// we want the one before that.
		if n > 0 {
			n--
		}
		ret = append(ret, float64(n))
	}
	return ret
}

type resAndChain struct {
	ResID int
	Chain string
}

// returns all the residue numbers in mol covered by indexes
func allResIDandChains(mol chem.Atomer, indexes []int) []*resAndChain {
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

// DistRank determines, for a reference set of coordinates and a set of residue names, the minimum distance between any atom from the reference
// and any atom from each residue with one of the names given (or the centroid of each residue, if a variadic "com" bool is given)
// returns a list with ID and distances, which satisfies the sort interface and has several other useful methods.
func DistRank(coord *v3.Matrix, mol chem.Atomer, refindexes []int, residues []string, options ...*Options) MolDistList {
	var o *Options
	if len(options) > 0 {
		o = options[0]
	} else {
		o = DefaultOptions()
	}

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
	cutoff := o.end
	cutoffplus := cutoff + 3
	//	chunk := NewTopology(0, 1)
	tmp := v3.Zeros(1)
	//	fmt.Println("Start looking!") ///////////////////
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		molname = at.MolName
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
			if o.com == true {
				var mass *mat.Dense
				topol := chem.NewTopology(0, 1)
				topol.SomeAtoms(mol, indexes)
				massp, err := topol.Masses()
				if err != nil {
					mass = mat.NewDense(len(massp), 1, massp)
				}
				c, err := chem.CenterOfMass(test, mass)
				if err != nil { //this really is very unlikely to fail. The ref matrix would have to be wrong. It could even deserve a panic.
					o.com = false //if it fails, we don't try again, for consistency. Any previous successful COM use will not comparable with numbers used from now on.
					log.Printf("gochem/solv/DistRank Couldn't obtain the COM/centroid: %s. Will work with all solvent atoms\n", err.Error())
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

// A structure for the distance from a residue to a particular point
type molDist struct {
	Distance float64
	MolID    int
}

func (M *molDist) str() string {
	return fmt.Sprintf("D: %4.3f ID: %d", M.Distance, M.MolID)
}

// A set of distances for different molecules to the same point
type MolDistList []*molDist

func (M MolDistList) Swap(i, j int) {
	M[i], M[j] = M[j], M[i]
}

// Less returns true if the distance of the element i
// to the pre-defined point is smallert than that of the element j,
// or false otherwise
func (M MolDistList) Less(i, j int) bool {
	return M[i].Distance < M[j].Distance
}
func (M MolDistList) Len() int {
	return len(M)
}

// Distance returns the distance from the element i of
// the slice to the pre-defined point
func (M MolDistList) Distance(i int) float64 {
	return M[i].Distance
}

// MolID resturns the MolID (the residue number, for a protein)
// of the i element of the slice
func (M MolDistList) MolID(i int) int {
	return M[i].MolID
}

// String produces a string representation of a set of distances
func (M MolDistList) String() string {
	retslice := make([]string, len(M))
	for i := range M {
		retslice[i] = M[i].str()
	}
	return strings.Join(retslice, "\n")
}

// Distances returns a slice with all the distances in the set
func (M MolDistList) Distances() []float64 {
	ret := make([]float64, len(M))
	for i := range M {
		ret[i] = M[i].Distance //This absolutely should never fail!
	}
	return ret
}

// MolIDs returns a slice with the molIDs in al the lists
func (M MolDistList) MolIDs() []int {
	ret := make([]int, len(M))
	for i := range M {
		ret[i] = M[i].MolID //This absolutely should never fail!
	}
	return ret
}

// Data returns a list with all the MolIDs and a list with all the distances
// in the set
func (M MolDistList) Data() ([]int, []float64) {
	ret := make([]int, len(M))
	ret2 := make([]float64, len(M))
	for i := range M {
		ret[i] = M[i].MolID //This absolutely should never fail!
		ret2[i] = M[i].Distance
	}
	return ret, ret2
}

// AtomIDs returns a list of all the atom IDs for all the residues in the list.
// Only a convenience function.
func (M MolDistList) AtomIDs(mol chem.Atomer) []int {
	molids := M.MolIDs()
	return chem.Molecules2Atoms(mol, molids, []string{})
}

// Merge aggregates the receiver and the arguments in the receiver
// it removes entries with repeated MolIDs
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

// MolShortestDistGiven two sets of coordinates, it obtains the distance between the two closest atoms
// in the 2 sets.
// This is probably not a very efficient way to do it.
// Note:This should probably be moved to the main gochem package, in geometric.go
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

//NOTE: These will be replaced when the generic funcions
//make it to Go's stdlib.

// isInInt returns true if test is in container, false otherwise.
func isInInt(container []int, test int) bool {
	if container == nil {
		return false
	}
	for _, i := range container {
		if test == i {
			return true
		}
	}
	return false
}

// Same as the previous, but with strings.
func isInString(container []string, test string) bool {
	if container == nil {
		return false
	}
	for _, i := range container {
		if test == i {
			return true
		}
	}
	return false
}
