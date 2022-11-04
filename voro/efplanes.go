package voro

import (
	"fmt"
	"sort"

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

func atomDistances(at int, c *v3.Matrix, cutoff float64, indexes ...[]int) dids {
	if cutoff < 0 {
		cutoff = defCutoff
	}
	//you can give this function a set of atoms to consider for the interactions
	in := make([]int, 0, 0)
	if len(indexes) == 0 || len(indexes[0]) == 0 {
		nvecs := c.NVecs()
		for i := 0; i < nvecs; i++ {
			in = append(in, i)
		}
	} else {
		in = indexes[0]
	}
	//	fmt.Println("in", len(in)) ////////////
	tmp := v3.Zeros(1)
	var test *v3.Matrix
	ref := c.VecView(at)
	vecs := c.NVecs()
	dists := make([]*did, 0, vecs-1)
	//This means many distances are calculated/stored twice.
	//I'll do things in this naive way for now
	for i := 0; i < vecs; i++ {
		if i == at || !isInInt(in, i) {
			continue
		}

		test = c.VecView(i)
		tmp.Sub(ref, test)
		d := tmp.Norm(2)
		if d < cutoff {
			dists = append(dists, &did{id: i, d: tmp.Norm(2)})
		}
		//else {
		//	fmt.Println("further away than cutoff! ", at, i, d, cutoff)
		//}

	}
	return dists
}

//if groups are given, each atom is only analized with respect to atoms in the group
//it is NOT present.
func distances(c *v3.Matrix, group1, group2 []int) []dids {
	vecs := c.NVecs()
	ret := make([]dids, 0, vecs)
	var only []int
	for i := 0; i < vecs; i++ {
		if group1 != nil && group2 != nil {
			if isInInt(group1, i) {
				only = group2
			} else if isInInt(group2, i) {
				only = group1
			}
		}
		//	fmt.Println("only", only) ////////////
		ret = append(ret, atomDistances(i, c, -1, only))
	}
	for _, v := range ret {
		sort.Sort(v)
	}
	return ret
}

//adds to the set of planes the plane from each atom to its (ith+1)th closest atom,
//if not repeated or blocked by already present planes.
func ithClosestPlanes(c *v3.Matrix, dists []dids, ith int, planes VPSlice) VPSlice {
	vecs := c.NVecs()
	blocked := 0
	repeated := 0
	total := 0
	excluded := 0
	for i := 0; i < vecs; i++ {
		//fmt.Println(len(dists[i]), ith) ///////////////////
		if len(dists[i]) <= ith {
			continue
		}
		//fmt.Println("found", ith, "atom to", i, "!") ///////////////////
		j := dists[i][ith].id
		ci := c.VecView(i)
		cj := c.VecView(j)
		p := PlaneBetweenAtoms(ci, cj, i, j)
		total++
		if planes.IsBlocked(p, c, i) { /////
			blocked++ /////////////
		} //////////////////////////////////
		if planes.IsRepeated(p) {
			repeated++
		}
		if planes.IsRepeated(p) || planes.IsBlocked(p, c, i) {
			excluded++
		}

		if !planes.IsBlocked(p, c, i) && !planes.IsRepeated(p) {
			planes = append(planes, p)
		}
	}
	fmt.Println(total, " total planes read, ", blocked, "planes blocked and ", repeated, "repeated at the ", "and", excluded, "planes excluded at the ", ith, "closest") ////////
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

func progressivePlanes(c *v3.Matrix, dists []dids) VPSlice {
	planess := make([]*VPlane, 0, 100)
	planes := VPSlice(planess)
	added := []int{1, 1, 1, 1} //, 1, 1, 1, 1}
	var add int
	//the following interates until no more planes are added
	for i := 0; !isInInt(added, 0); i++ {
		println("vueeeltaaa") ///////
		prev := len(planes)
		planes = ithClosestPlanes(c, dists, i, planes)
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

func ContactPlanes(c *v3.Matrix, group1, group2 []int) VPSlice {
	dists := distances(c, group1, group2)
	atoms := 0
	avdist := 0
	for _, v := range dists {
		avdist += len(v)
		atoms++
	}
	fmt.Println("average distances", avdist/atoms) ////////////
	//fmt.Println(dists) /////////////////////
	return progressivePlanes(c, dists)
}
