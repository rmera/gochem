/*
 *   areas.go, part of gochem.
 *   Obtain areas for the faces of voronoi polihedra
 *
 *
 *
 */
package voro

import (
	"fmt"
	"math"

	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

type edge struct {
	plances  []*VPlane //Pointer to at least 2 of the planes that intersect in this edge
	vertices []*vertix
}

type vertix struct {
	planes   []*VPlane //pointer to the planes intersecting on this vertix
	loc      *v3.Matrix
	isVertix bool
}

func findVertix(p1, p2, p3 *VPlane, Adata, Ainvdata, Cdata []float64, pars [12]float64) *vertix {
	if len(Adata) < 9 || len(Ainvdata) < 9 || len(Cdata) < 3 {
		s := fmt.Sprintf("len Adata: %d, len Ainvdata: %d, len Cdata: %d\n", len(Adata), len(Ainvdata), len(Cdata))
		panic("Was given wrong temporary slices: " + s)
	}
	tempdense := mat.NewDense(1, 3, Cdata)
	temp := v3.Dense2Matrix(tempdense)
	checkplanes := []*VPlane{p1, p2, p3}
	//we first check that the solution exists
	for i, v := range checkplanes {
		for _, w := range checkplanes[i+1:] {
			temp.Cross(v.Normal, w.Normal)
			if temp.Norm(2) == 0 {
				return nil
			}
		}
	}
	p1.Parametric(pars[:4])
	p2.Parametric(pars[4:8])
	p3.Parametric(pars[8:12])
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			Adata[3*i+j] = pars[3*i+j+i] //we always skip the i*3 element in pars (goes to the C  matrix)
			//hence the additional "i" that is added to that index
		}
	}
	A := mat.NewDense(3, 3, Adata)
	Cdata[0] = pars[3]
	Cdata[1] = pars[7]
	Cdata[2] = pars[11]
	C := mat.NewDense(3, 1, Cdata)
	Ainv := mat.NewDense(3, 3, Ainvdata)
	err := Ainv.Inverse(A)
	if err != nil {
		panic("gonum/mat/Dense.inverse: " + err.Error())
	}
	res := mat.NewDense(3, 1, make([]float64, 3))
	res.Mul(Ainv, C)
	resdata := res.RawMatrix().Data
	resT := mat.NewDense(1, 3, resdata)
	return &vertix{planes: []*VPlane{p1, p2, p3}, loc: v3.Dense2Matrix(resT), isVertix: true}
}

//isInInt is a helper for the RamaList function,
//returns true if test is in container, false otherwise.
//At some point I'll upgrade to generics :-D
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

func closest(temp *v3.Matrix, target int, excluded []int, verts []*vertix) int {
	closest := 0
	dist := 0.0
	dtemp := 0.0
	for i, v := range verts {
		if i == target || isInInt(excluded, i) {
			continue
		}
		temp.Sub(verts[target].loc, v.loc)
		dtemp = temp.Norm(2)
		if dtemp < dist || i == 0 {
			closest = i
			dist = dtemp
		}
	}
	return closest

}

func herontriangleArea(p1, p2, o, tmp *v3.Matrix, sides []*v3.Matrix) float64 {
	var a, b, c, s float64
	tmp.Sub(sides[0], sides[1])
	a = tmp.Norm(2)
	tmp.Sub(sides[1], sides[2])
	b = tmp.Norm(2)
	tmp.Sub(sides[2], sides[0])
	c += tmp.Norm(2)
	s = (a + b + c) / 2.0
	return math.Sqrt(s * (s - a) * (s - b) * (s - c))
}

//it is assumed that only the vertices belonging to the face in question are included in verts
func areaFromVertices(atomCoord *v3.Matrix, verts []*vertix, tmp *v3.Matrix) float64 {
	ordered := make([]int, 0, 0)
	temp := v3.Zeros(1)
	fmt.Println("largos ql", len(verts), len(ordered))
	for i := 0; len(ordered) < len(verts); i = closest(temp, i, ordered, verts) {
		fmt.Println("vertice", i) ///////////////
		ordered = append(ordered, i)
	}
	area := 0.0
	sides := []*v3.Matrix{v3.Zeros(1), v3.Zeros(1), v3.Zeros(1)}
	last := len(ordered) - 1
	for _, v := range ordered[:last] { //all but the last triangle
		area += herontriangleArea(verts[v].loc, verts[v+1].loc, atomCoord, tmp, sides)
	}
	//the last triangle
	area += herontriangleArea(verts[ordered[last]].loc, verts[ordered[0]].loc, atomCoord, tmp, sides)
	return area

}

//PairContactArea takes the vertices for the voronoi polihedra of an atom at1, plus the index for that
//atom and that of a second atom at2. It returns the area for the face of the voronoi polihedron separating
//the 2 atoms.
func PairContactArea(at1, at2 int, coords *v3.Matrix, planes VPSlice) float64 {
	verts := verticesInAtomPair(planes, at1, at2, coords.VecView(at1))
	fmt.Println(verts) /////////////
	if len(verts) == 0 {
		println("so no vertices, it seems")
		return 0 //It might be that in some cases, the area between 2 atoms is "unbound" so you can't estimate it like this.
	}
	//we first need to select the vertices that are shared with at2
	tmp := v3.Zeros(1)
	return areaFromVertices(coords.VecView(at1), verts, tmp)
}

func sameAtoms(a, b []int) bool {
	if a[0] == b[0] && a[1] == b[1] {
		return true
	}
	return false
}

//Looks for al the vertices for the voronoi polihedra which are in the plane separating atoms
//atom and atom2, i.e. the limits of the face of the polihedron that separates those 2
//atoms.
func verticesInAtomPair(planes VPSlice, atom, atom2 int, atomCoord *v3.Matrix) []*vertix {
	ref := planes.PairPlane(atom, atom2)
	var pars [12]float64
	Adata := make([]float64, 9)
	Ainvdata := make([]float64, 9)
	Cdata := make([]float64, 3)
	planesat := filterPerAtom(planes, atom)
	//fmt.Println("planesat", len(planesat)) //////////
	ret := make([]*vertix, 0, 5)
	//all possible triads of planes
	for i, v := range planesat {
		for _, w := range planesat[i+1:] {
			if sameAtoms(v.Atoms, ref.Atoms) || sameAtoms(w.Atoms, ref.Atoms) {
				continue
			}
			vertix := findVertix(v, w, ref, Adata, Ainvdata, Cdata, pars)
			if vertix != nil {
				ret = append(ret, vertix)
			}
		}
	}
	//We now have all vertices between planes, but only some of them are actual vertices of the polihedron.
	//In many cases the vertix will be "blocked" by another of the planes.
	//We
	od := v3.Zeros(1)
	var dvert, dp float64
	total := len(ret)
	//println(len(ret), "total vertices") /////////////////
	for _, v := range ret {
		od.Sub(v.loc, atomCoord)
		for _, p := range planesat {
			if p == v.planes[0] || p == v.planes[1] || p == v.planes[2] {
				continue //as the planes are all pointers, this should check that they point to the same address, i.e.
				//are the same. There is no point in checking the planes that form the vertices.
			}
			od.Sub(v.loc, atomCoord)
			dp = p.DistanceInterVector2(atomCoord, v.loc)
			if dp < 0 {
				continue //it means the plane never intesects the vector
			}
			dvert = od.Norm(2)
			if dp < dvert {
				///	fmt.Println("falsalaweaaaaaaaa", dp, dvert) //////
				v.isVertix = false
				total--
				break
			}
		}

	}
	retdef := make([]*vertix, 0, total)
	for _, v := range ret {
		if v.isVertix {

			retdef = append(retdef, v)
		}
	}
	fmt.Println("vertices finales!", len(retdef), total) ///////////

	return retdef

}

func filterPerAtom(planes VPSlice, atom int) VPSlice {
	ret := make([]*VPlane, 0, 5)
	for _, v := range planes {
		if v.Atoms[0] == atom || v.Atoms[1] == atom {
			ret = append(ret, v)
		}
	}
	return ret
}
