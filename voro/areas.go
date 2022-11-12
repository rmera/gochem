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
	"log"
	"math"

	chem "github.com/rmera/gochem"
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
	dist := 99999999999.0
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

func herontriangleArea(p1, p2, o, tmp *v3.Matrix) float64 {
	var a, b, c, s float64

	tmp.Sub(o, p1)
	a = tmp.Norm(2)
	tmp.Sub(p1, p2)
	b = tmp.Norm(2)
	tmp.Sub(p2, o)
	c = tmp.Norm(2)
	s = (a + b + c) / 2.0
	return math.Sqrt(s * (s - a) * (s - b) * (s - c))
}

func areaAndVolumeFromVerticesAndPlane(atomCoord *v3.Matrix, verts []*vertix, plane *VPlane, tmp *v3.Matrix) (area float64, volume float64) {

	ordered := make([]int, 0, 0)
	temp := v3.Zeros(1)
	center := v3.Zeros(1)
	for _, v := range verts {
		center.Add(center, v.loc)
	}
	center.Scale(1/float64(len(verts)), center)
	//fmt.Println("largos ql", len(verts), len(ordered))
	for i := 0; len(ordered) < len(verts); i = closest(temp, i, ordered, verts) {
		//	fmt.Println("vertice", i) ///////////////
		ordered = append(ordered, i)
	}
	last := len(ordered) - 1
	for i, v := range ordered[:last] { //all but the last triangle
		a := herontriangleArea(verts[v].loc, verts[ordered[i+1]].loc, center, tmp)
		h := math.Abs(plane.DistanceInterVector(atomCoord, plane.Normal))
		area += a
		volume += (1 / 3.0) * a * h

	}
	//the last triangle
	area += herontriangleArea(verts[ordered[last]].loc, verts[ordered[0]].loc, center, tmp)
	return area, volume

}

//it is assumed that only the vertices belonging to the face in question are included in verts
func areaFromVertices(atomCoord *v3.Matrix, verts []*vertix, tmp *v3.Matrix) float64 {
	ordered := make([]int, 0, 0)
	temp := v3.Zeros(1)
	center := v3.Zeros(1)
	for _, v := range verts {
		center.Add(center, v.loc)
	}
	center.Scale(1/float64(len(verts)), center)
	//fmt.Println("largos ql", len(verts), len(ordered))
	for i := 0; len(ordered) < len(verts); i = closest(temp, i, ordered, verts) {
		//	fmt.Println("vertice", i) ///////////////
		ordered = append(ordered, i)
	}
	area := 0.0
	last := len(ordered) - 1
	for i, v := range ordered[:last] { //all but the last triangle
		//	fmt.Println("center, vertices", verts[v].loc, verts[ordered[i+1]].loc, center) ////
		area += herontriangleArea(verts[v].loc, verts[ordered[i+1]].loc, center, tmp)
		//	fmt.Println(verts[v], verts[ordered[i+1]], v, ordered[i+1], "vertices ql", area, "area ql") /////////////////////
	}
	//the last triangle
	area += herontriangleArea(verts[ordered[last]].loc, verts[ordered[0]].loc, center, tmp)
	return area

}

//PairContactAreaAndVolume takes the vertices for the voronoi polihedra of an atom at1, plus the index for that
//atom and that of a second atom at2. It returns the area for the face of the voronoi polihedron separating the 2 atoms,
func PairContactArea(at1, at2 int, coords *v3.Matrix, mol chem.Atomer, planes VPSlice) float64 {
	ret := PairContactAreaAndVolume(at1, at2, coords, mol, planes, false)
	return ret[0]
}

//PairContactAreaAndVolume takes the vertices for the voronoi polihedra of an atom at1, plus the index for that
//atom and that of a second atom at2. It returns an array of 2 float64 where the first element is the area for the face of the voronoi polihedron separating the 2 atoms, and the second is 0. If volume is given (only the first element of the slice is considered) and it is true, the volume of the polyhedron centerd in at1 associated to at2 is also obtained and returned as the second element of the slice.
func PairContactAreaAndVolume(at1, at2 int, coords *v3.Matrix, mol chem.Atomer, planes VPSlice, notvolume ...bool) []float64 {
	verts := verticesInAtomPair(planes, at1, at2, coords.VecView(at1))
	if len(verts) < 3 {
		log.Println("The interface between toms", at1, at2, "seem to be unbounded, or they could be sharing only an edge")
		return []float64{0, 0} //It might be that in some cases, the area between 2 atoms is "unbound" so you can't estimate it like this.
	}
	//	vertstr := ""
	//	for i, v := range verts {
	//		vertstr += fmt.Sprintf("n %d loc: %v |", i, v.loc)
	//	}
	//	println(vertstr) //////////////////
	//we first need to select the vertices that are shared with at2
	tmp := v3.Zeros(1)
	if len(notvolume) > 0 && notvolume[0] {
		area := areaFromVertices(coords.VecView(at1), verts, tmp)
		return []float64{area, 0}

	}

	plane := planes.PairPlane(at1, at2)
	area, volume := areaAndVolumeFromVerticesAndPlane(coords.VecView(at1), verts, plane, tmp)
	vdw := mol.Atom(at1).Vdw
	vdw2 := mol.Atom(at2).Vdw
	if vdw2 < vdw {
		vdw = vdw2
	}
	//I am assuming that the contact area should never be greater than the 'disk' are of the hemisphere
	//of the smaller atom. This probably means at1 is unbounded
	if area > math.Pi*math.Pow(vdw, 2) {
		log.Printf("Too large area of %3.5f A for contact between atoms %d and %d. Area of the 'disk' in the atom's hemispheres are %3.5f and %3.5f A^2, the smallest of which should be the maximum possible contact area. Likely unbounded polyhedron for atom %d, the area and volume will be set to 0.0\n", area, at1, at2, math.Pi*math.Pow(mol.Atom(at1).Vdw, 2), math.Pi*math.Pow(mol.Atom(at2).Vdw, 2), at1)
		return []float64{0, 0}
	}

	return []float64{area, volume}

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
	planesat := planes.AtomPlanes(atom) //filterPerAtom(planes, atom)
	//planesat = append(planesat, filterPerAtom(planes, atom2)...) //////////////////////////////////
	//	println("planesat", len(planesat)) /////////
	if len(planesat) < 4 {
		return nil
	}
	ret := make([]*vertix, 0, 5)
	//all possible triads of planes
	deflongvertex := false
	for i, v := range planesat {
		for _, w := range planesat[i+1:] {
			if sameAtoms(v.Atoms, ref.Atoms) || sameAtoms(w.Atoms, ref.Atoms) {
				//fmt.Println("skipped", i, v.Atoms, w.Atoms, ref.Atoms) ///////
				continue
			}
			vertix := findVertex(v, w, ref, Adata, Ainvdata, Cdata, pars)
			if vertix != nil {
				ret = append(ret, vertix)
			}
		}
	}
	longvertex := false
	//We now have all vertices between planes, but only some of them are actual vertices of the polihedron.
	//In many cases the vertix will be "blocked" by another of the planes.
	//We
	od := v3.Zeros(1)
	var dvert, dp float64
	total := len(ret)
	for _, v := range ret {
		od.Sub(v.loc, atomCoord)
		dvert = od.Norm(2)
		if dvert > 50 { /////////////
			longvertex = true
			//			fmt.Println("vertice largo", dvert, "planesat:", len(planesat))
		}
		for _, p := range planesat {
			if p == v.planes[0] || p == v.planes[1] || p == v.planes[2] {
				//println("pasa?") //////////////////////////
				continue //as the planes are all pointers, this should check that they point to the same address, i.e.
				//are the same. There is no point in checking the planes that form the vertices.
			}
			dp = p.DistanceInterVector(atomCoord, v.loc)
			//			if dvert > 50 {
			//				println("a plane's distance", dp) //////////
			//			}

			//println("distancia reql", dp) ///////////////////
			if dp < 0 {
				//				if dvert > 50 {
				//					println("PassedDue dp<0", dp) //////////
				//				}
				continue //it means the plane never intesects the vector
			}
			if dp < dvert {
				longvertex = false
				//			if dvert > 50 {
				//				fmt.Println("nos salvamos!", dvert, dp, p.Atoms) //////////

				//			}

				//println(dp, dvert) ////////////////////////////
				v.isVertix = false
				total--
				break
			}
		}
		if longvertex {
			deflongvertex = true
		}

	}
	retdef := make([]*vertix, 0, total)
	for _, v := range ret {
		if v.isVertix {

			retdef = append(retdef, v)
		}
	}
	//debugging stuff, draw the long vertices in an xyz file
	_ = deflongvertex
	if len(retdef) > 3 {
		cvert := v3.Zeros(len(retdef))
		ats := make([]*chem.Atom, len(retdef))
		filename := fmt.Sprintf("vert_%d-%d.xyz", atom, atom2)
		for i, v := range retdef {
			ats[i] = &chem.Atom{Symbol: "I"}
			cvert.Set(i, 0, v.loc.At(0, 0))
			cvert.Set(i, 1, v.loc.At(0, 1))
			cvert.Set(i, 2, v.loc.At(0, 2))
		}
		mol2 := chem.NewTopology(0, 1, ats)
		chem.XYZFileWrite(filename, cvert, mol2)

	}
	//end debugging stuff
	//	fmt.Println("original and cleaned", len(ret), len(retdef))
	return retdef

}

//all the "data" and "pars" arguments are things to be used for memory by findvertex. The data theyc ontain will be errased.
func findVertex(p1, p2, p3 *VPlane, Adata, Ainvdata, Cdata []float64, pars [12]float64) *vertix {
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
		return nil //there is jut no solution
		//panic("gonum/mat/Dense.inverse: " + err.Error())
	}
	res := mat.NewDense(3, 1, make([]float64, 3))
	res.Mul(Ainv, C)
	resdata := res.RawMatrix().Data
	resT := mat.NewDense(1, 3, resdata)
	return &vertix{planes: []*VPlane{p1, p2, p3}, loc: v3.Dense2Matrix(resT), isVertix: true}
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
