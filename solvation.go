package chem

import (
	"fmt"
	"runtime"
	"strings"

	//	"sort"
	//	"strconv"
	v3 "github.com/rmera/gochem/v3"
)

func SerialMolRDF(traj Traj, mol *Molecule, refindexes []int, residues []string, step, end float64, frameskip int) ([]float64, error) {
	var ret []float64
	coords := v3.Zeros(mol.Len())
	framesread := 0
reading:
	for i := 0; ; i++ {
		err := traj.Next(coords)
		if err != nil {
			switch err := err.(type) {
			default:
				return nil, err
			case LastFrameError:
				break reading
			}
		}
		rdf := FrameUMolCRDF(coords, mol, refindexes, residues, step, end)
		if ret == nil {
			ret = make([]float64, len(rdf))
		}
		for j, _ := range ret {
			ret[j] = ret[j] + rdf[j]
		}
		framesread += 1

	}
	//yeeeeah, I'm not super sure about this.
	avtotalsolv := ret[len(ret)-1] / float64(framesread)
	for i, _ := range ret {
		if i > 0 {
			ret[i] = ret[i] - ret[i-1]
		}
		ret[i] = ret[i] / avtotalsolv
		ret[i] = ret[i] / float64(framesread)

	}

	return ret, nil
}

//Calculates the RDF for a trajectory given the indexes of the solute atoms, the solvent molecule name, the step for the "layers" and the cutoff.
func MolRDF(traj ConcTraj, mol *Molecule, refindexes []int, residues []string, step, end float64, frameskip int) ([]float64, error) {
	var ret []float64
	simulframes := runtime.NumCPU()
	runtime.GOMAXPROCS(simulframes) //shouldn't be needed anymore, I think.
	frames := make([]*v3.Matrix, simulframes, simulframes)
	framesread := 0
	for i, _ := range frames {
		frames[i] = v3.Zeros(traj.Len())
	}
	results := make([][]chan []float64, 0, 0)
	for i := 0; ; i++ {
		results = append(results, make([]chan []float64, 0, len(frames)))
		coordchans, err := traj.NextConc(frames)
		if err != nil {
			if err, ok := err.(LastFrameError); ok {
				if coordchans == nil {
					break
				}
			} else {
				return nil, err //Must decorate this error
			}
		}
		for key, channel := range coordchans {
			results[len(results)-1] = append(results[len(results)-1], make(chan []float64))
			go unitRDF(channel, results[len(results)-1][key], mol, refindexes, residues, step, end)
		}
		res := len(results) - 1
		for _, k := range results[res] {
			if k == nil {
				continue //this means we skipped the frame. Right now, it should never happen.
			}
			framesread++
			rettemp := <-k
			if ret == nil {
				ret = make([]float64, len(rettemp))
			}
			for i, _ := range ret {
				ret[i] = ret[i] + rettemp[i]
			}
		}
	}
	//yeeeeah, I'm not super sure about this.
	avtotalsolv := ret[len(ret)-1] / float64(framesread)
	for i, _ := range ret {
		if i > 0 {
			ret[i] = ret[i] - ret[i-1]
		}
		ret[i] = ret[i] / avtotalsolv
		ret[i] = ret[i] / float64(framesread)

	}

	return ret, nil
}

//The worker function for the RDF
func unitRDF(channelin chan *v3.Matrix, channelout chan []float64, mol *Molecule, refindexes []int, residues []string, step, end float64) {
	if channelin != nil {
		temp := <-channelin
		rdf := FrameUMolCRDF(temp, mol, refindexes, residues, step, end)
		channelout <- rdf
	} else {
		channelout <- nil
	}
	return
}

//obtains the Unnormalized "Cummulative Molecular RDF" for one solvated structure. The RDF would be these values averaged over several structures.
func FrameUMolCRDF(coord *v3.Matrix, mol *Molecule, refindexes []int, residues []string, step, end float64) []float64 {
	if step <= 0 {
		step = 0.5
	}
	if end <= 0 {
		end = 15
	}
	totalsteps := int(end/step) + 1
	ret := make([]float64, 0, totalsteps)
	// Here we get the RDF by counting the molecules at X distance or less from the solute, and subtracting the result from the previous callculation.
	//This means that we do more calculation than needed, as every time keep including the solvent that was in previous layers.
	//	acclen := 0.0
	for i := 0.0; i <= end; i = i + step {
		result := DistRank(coord, mol, refindexes, residues, i)
		l := float64(result.Len())
		//		acclen = l //keeps track of the last calculated. To be used for later normalization
		//		if len(ret) >= 1 {
		//			l = l - ret[len(ret)-1]
		//		}
		ret = append(ret, l)
	}
	return ret
}

//obtains the unnormalized "Cummulative Molecular RDF" for one solvated structure. The RDF would be these values averaged over several structures.
//works concurrently
func ConcFrameUMolCRDF(coord *v3.Matrix, mol *Molecule, refindexes []int, residues []string, step, end float64) []float64 {
	if step <= 0 {
		step = 0.5
	}
	if end <= 0 {
		end = 15
	}
	totalsteps := int(end/step) + 1
	ret := make([]float64, totalsteps)
	results := make([][]chan float64, 0, totalsteps)

	simulchunks := runtime.NumCPU()
	runtime.GOMAXPROCS(simulchunks) //shouldn't be needed anymore, I think.
	for i := 0; i < totalsteps; i = i + simulchunks {
		results = append(results, make([]chan float64, 0, simulchunks))
		reslen := len(results) - 1
		for j := 0; j < simulchunks || i+j <= totalsteps; j++ {
			results[reslen] = append(results[reslen], make(chan float64))
			go func(reschan chan float64, cutoff float64) {

				res := DistRank(coord, mol, refindexes, residues, cutoff)
				reschan <- float64(res.Len())

			}(results[reslen][i], float64(i+j)*step)
		}
		//now get the results from the workers
		for index, val := range results[reslen] {
			ret[i+index] = <-val
		}

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

func DistRank(coord *v3.Matrix, mol *Molecule, refindexes []int, residues []string, cutoff float64) MolDistList {
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
			for j := i + 1; j < i+80; j++ {
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
			distance = MolShortestDist(test, ref)
			if distance <= cutoff {
				ranks = append(ranks, &molDist{Distance: distance, MolID: id})
				molid_skip = id
				chain_skip = chain
			}
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

//MolShortestDistGiven two sets of coordinates, it obtains the shortest distance from any 2 points
//in the set. This is probably not a very efficient way to do it.
func MolShortestDist(test, ref *v3.Matrix) float64 {
	temp := v3.Zeros(1)
	var d1, dclosest float64
	var vt1, vr1 *v3.Matrix // vtclosest,vr1, vrclosest *v3.Matrix
	dclosest = 100000       //we initialize this to some crazy large value so it's immediately replaced with the first calculated distance.
	// dist(ref.VecView(0), test.VecView(0), temp)
	for i := 0; i < test.NVecs(); i++ { //We do on not-needed comparison
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
