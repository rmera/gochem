package voro

import (
	"fmt"
	"sort"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

type did struct {
	d  float64
	id int
}

type dids []*did

func (d dids) Less(i, j int) bool {
	return d[i].d < d[j].d
}
func (d dids) Len() int {
	return len(d)
}
func (d dids) Swap(i, j int) {
	d[i], d[j] = d[j], d[i]
}

func (d dids) String() string {
	ret := ""
	for _, v := range d {
		ret += fmt.Sprintf(" | id: %d ; d: %5.3f", v.id, v.d)
	}
	return ret
}

func atomVdwDistances(at int, c *v3.Matrix, mol chem.Atomer, noH bool, allowed []int, scan ...*AngleScan) []*did {
	if len(scan) <= 0 {
		scan = append(scan, DefaultAngleScan()) //1.4 is the vdW radius of water

	}
	var vdwsum, cutoff float64
	tmp := v3.Zeros(1)
	var test *v3.Matrix
	ref := c.VecView(at)
	vecs := c.NVecs()
	dists := make([]*did, 0, vecs-1)
	//This means many distances are calculated/stored twice.
	//I'll do things in this naive way for now
	for i := 0; i < vecs; i++ {
		if i == at || !isInInt(allowed, i) {
			continue
		}
		if noH && mol.Atom(i).Symbol == "H" {
			continue
		}

		vdwsum = mol.Atom(at).Vdw + mol.Atom(i).Vdw
		cutoff = vdwsum*scan[0].VdwFactor + scan[0].Offset

		test = c.VecView(i)
		tmp.Sub(ref, test)
		d := tmp.Norm(2)
		if d < cutoff {
			dists = append(dists, &did{id: i, d: d})
		}
		//else {
		//	fmt.Println("further away than cutoff! ", at, i, d, cutoff)
		//}

	}
	return dists

}

func vdwdistances(c *v3.Matrix, mol chem.Atomer, noH bool, alo []int) []dids {
	vecs := c.NVecs()
	ret := make([]dids, 0, vecs)
	for i := 0; i < vecs; i++ {
		//	fmt.Println("only", only) ////////////
		if !isInInt(alo, i) {
			ret = append(ret, nil)
			continue
		}
		if noH && mol.Atom(i).Symbol == "H" {
			ret = append(ret, nil)
			continue
		}
		ret = append(ret, atomVdwDistances(i, c, mol, noH, alo))
	}
	for _, v := range ret {
		sort.Sort(v)
	}
	return ret
}

//adds to the set of planes the plane from each atom to its (ith+1)th closest atom,
//if not repeated or blocked by already present planes.
func kthClosestPlanes(c *v3.Matrix, mol chem.Atomer, dists []dids, kth int, noH bool, allowed []int, planes VPSlice) VPSlice {
	vecs := c.NVecs()
	total := 0
	offset := 0
	for i := 0; i < vecs; i++ {
		if !isInInt(allowed, i) {
			continue
		}
		if noH && mol.Atom(i).Symbol == "H" {
			continue
		}
		//fmt.Println(len(dists[i]), ith) ///////////////////
		if len(dists[i]) <= kth {
			//			println("not enough distances!!1") ////////
			continue
		}
		//fmt.Println("found", ith, "atom to", i, "!") ///////////////////
		j := dists[i][kth].id
		ci := c.VecView(i)
		cj := c.VecView(j)
		p := PlaneBetweenAtoms(ci, cj, i, j)
		total++
		//		if planes.IsBlocked(p, c, i) { /////
		//			blocked++ /////////////
		//		} //////////////////////////////////
		//		if planes.IsRepeated(p) {
		//			repeated++
		//		}
		//		if planes.IsRepeated(p) || planes.IsBlocked(p, c, i) {
		//			excluded++
		//		}

		if !planes.IsRepeated(p) && !planes.IsBlocked(p, c, i-offset) {
			planes = append(planes, p)
		}
	}
	fmt.Println(total, " total planes read") ////////
	//now we prune the repeated planes
	newplanes := make([]*VPlane, 0, len(planes)/2)
	np := VPSlice(newplanes)
	for _, v := range planes {
		if !np.IsRepeated(v) {
			np = append(np, v)
		}
	}

	return planes
}

func progressivePlanes(c *v3.Matrix, mol chem.Atomer, dists []dids, noH bool, allowed []int) VPSlice {
	planess := make([]*VPlane, 0, 100)
	planes := VPSlice(planess)
	added := []int{1, 1} //, 1, 1, 1, 1}
	var add int
	//the following interates until no more planes are added
	for i := 0; !isInInt(added, 0); i++ {
		println("vueeeltaaa", i) ///////
		prev := len(planes)
		planes = kthClosestPlanes(c, mol, dists, i, noH, allowed, planes)
		post := len(planes)
		add = post - prev
		for j := 0; j < len(added)-1; j++ {
			added[j] = added[j+1]
		}
		added[len(added)-1] = add
		fmt.Println("planes post, prev, added:", post, prev, added) ///////////
	}
	return planes
}

func ContactPlanes(c *v3.Matrix, mol chem.Atomer, noH bool, allowed ...[]int) VPSlice {
	var alo []int
	vecs := c.NVecs()
	if len(allowed) > 0 && len(allowed[0]) > 0 {
		alo = allowed[0]
	} else {
		alo = make([]int, 0, vecs)
		for i := 0; i < vecs; i++ {
			alo = append(alo, i)
		}

	}

	dists := vdwdistances(c, mol, noH, alo)
	atoms := 0
	avdist := 0
	for _, v := range dists {
		avdist += len(v)
		atoms++
	}
	fmt.Println("average distances", avdist/atoms) ////////////
	//fmt.Println(dists) /////////////////////
	return progressivePlanes(c, mol, dists, noH, alo)
}
