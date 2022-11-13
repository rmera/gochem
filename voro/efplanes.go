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

func atomVdwDistances(at int, c *v3.Matrix, mol chem.Atomer, scan *Options) []*did {
	var vdwsum, cutoff float64
	tmp := v3.Zeros(1)
	var test *v3.Matrix
	ref := c.VecView(at)
	vecs := c.NVecs()
	dists := make([]*did, 0, vecs-1)
	//This means many distances are calculated/stored twice.
	//I'll do things in this naive way for now
	for i := 0; i < vecs; i++ {
		if i == at || !isInInt(scan.Subset, i) {
			continue
		}
		if scan.NoH && mol.Atom(i).Symbol == "H" {
			continue
		}
		vdwsum = mol.Atom(at).Vdw + mol.Atom(i).Vdw + scan.vdwSum
		cutoff = vdwsum*scan.VdwFactor + scan.Offset

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

func vdwdistances(c *v3.Matrix, mol chem.Atomer, scan *Options) []dids {
	vecs := c.NVecs()
	ret := make([]dids, 0, vecs)
	alo := scan.Subset
	for i := 0; i < vecs; i++ {
		//	fmt.Println("only", only) ////////////
		if !isInInt(alo, i) {
			ret = append(ret, nil)
			continue
		}
		if scan.NoH && mol.Atom(i).Symbol == "H" {
			ret = append(ret, nil)
			continue
		}
		tmp := atomVdwDistances(i, c, mol, scan)
		sort.Sort(dids(tmp))
		ret = append(ret, tmp)
	}
	return ret
}

//adds to the set of planes the plane from each atom to its (ith+1)th closest atom,
//if not repeated or blocked by already present planes.
func kthClosestPlanes(c *v3.Matrix, mol chem.Atomer, dists []dids, kth int, planes VPSlice, options ...*Options) VPSlice {
	vecs := c.NVecs()
	total := 0
	//repeated := 0
	//blocked := 0
	for i := 0; i < vecs; i++ {
		if len(dists[i]) <= kth {
			continue
		}
		j := dists[i][kth].id
		ci := c.VecView(i)
		cj := c.VecView(j)
		//Unless the options request it, the plane is no longer in the middle of both atoms, but in the weighted average
		//between them, where the weights are the vdW radii.
		vdw1 := mol.Atom(i).Vdw
		vdw2 := mol.Atom(j).Vdw
		if len(options) >= 1 && options[0].EquidistantPlanes {
			vdw1 = 1.0
			vdw2 = 1.0
		}
		p := PlaneBetweenAtoms(ci, cj, i, j, vdw1, vdw2)
		total++
		if !planes.IsRepeated(p) && !planes.IsBlocked(p, c, i) {
			planes = append(planes, p)
		}
	}
	//	fmt.Println(total, " total planes read. ", repeated, "found repeated", blocked, "found blocked") ////////
	return planes
}

func progressivePlanes(c *v3.Matrix, mol chem.Atomer, dists []dids, options ...*Options) VPSlice {
	planess := make([]*VPlane, 0, 100)
	planes := VPSlice(planess)
	added := []int{1, 1}
	var add int
	//the following interates until no more planes are added
	for i := 0; !isInInt(added, 0); i++ {
		//println("vueeeltaaa", i) ///////
		prev := len(planes)
		planes = kthClosestPlanes(c, mol, dists, i, planes, options...)
		post := len(planes)
		add = post - prev
		for j := 0; j < len(added)-1; j++ {
			added[j] = added[j+1]
		}
		added[len(added)-1] = add
		//fmt.Println("planes post, prev, added:", post, prev, added) ///////////
	}
	return planes
}

func ContactPlanes(c *v3.Matrix, mol chem.Atomer, options ...*Options) VPSlice {
	var scan *Options
	if len(options) == 0 || options[0] == nil {
		scan = DefaultOptions()
	} else {
		scan = options[0]
	}
	if scan.NoH {
		scan.vdwSum = 2
	}
	if scan.WaterContacts {
		scan.Offset = defWaterOffset
	}
	var alo []int
	vecs := c.NVecs()
	if len(scan.Subset) == 0 {
		scan.Subset = make([]int, 0, vecs)
		for i := 0; i < vecs; i++ {
			scan.Subset = append(alo, i)
		}

	}
	dists := vdwdistances(c, mol, scan)
	atoms := 0
	avdist := 0
	for _, v := range dists {
		avdist += len(v)
		atoms++
	}
	//fmt.Println("average distances", avdist/atoms) ////////////
	//fmt.Println(dists) /////////////////////
	return progressivePlanes(c, mol, dists, options...)
}
