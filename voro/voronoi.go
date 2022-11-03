/*
 * voronoi.go, part of gochem.
 *
 *
 * Copyright 2021 Raul Mera <rauldotmeraatusachdotcl>
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
 * goChem is currently developed at the Universidad de Santiago de Chile
 * (USACH)
 *
 */

/**
voro offers Voronoi polyhedra-based functionality for goChem, allowing to determine whether 2 atoms
are in contact.

Unfortunately, due to the author's lack of knowledge on comptuational geometry, and the lack of implementations
of 3D-Voronoi polyhedra construction in Go, this library is extremely naive, un-sophisticated and
heavily brute-force, in addition to incomplete. It uses a series of rather dirty tricks invoking physical and
chemical knowledge to alliviate the performance problems caused by the previous, but they are still here.
Also, the heavy numerical character of voro is sure to cause inaccuracies when compared to a proper,
analytical implementation. Still, it does work in my tests.

**/

package voro

import (
	"fmt"
	"log"
	"os"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/stat/combin"
)

const (
	defOffset    float64 = 1.4
	defVdwFactor float64 = 1.2
	defMaxAngle  float64 = 87
	defAngleStep float64 = 5
	defCutoff    float64 = 6 //rather permissive
)

//AngleScan contains options to perform angle scans to see if there is an angle in which 2 atoms
//are in direct contact (i.e. if part of the plane bisecting them is part of the Voronoi polihedra for the system).
type AngleScan struct {
	VdwFactor float64
	Offset    float64
	Angles    []float64 //last angle, step between angles, in degrees. A full scan would be 0 to 90
	Cutoff    float64   //distance cutoff, if the vectors at the given angle are farther than this, the angle is ignored.
	Test      bool      //for debugging
}

//DefaultAngleScan returns the default setting for an AngleScan
func DefaultAngleScan() *AngleScan {
	return &AngleScan{Offset: defOffset, VdwFactor: defVdwFactor, Angles: []float64{defMaxAngle, defAngleStep}, Cutoff: defCutoff}
}

//This is a naive, unoptimal, simple and incomplete implementation
//of 3D Voronoi polihedra.

//The plane bisecting 2 atoms
//With all the info we might want from it.
type VPlane struct {
	Atoms      []int   //indexes of the 2 atoms.
	Distance   float64 //The plane is equidistant from both atoms, so there is only one distance
	Normal     *v3.Matrix
	Point      *v3.Matrix //A point in the plane.
	NotContact bool       //true if this plane is a confirmed non-contact
	Contact    bool       //true if this plane is a confirmed contact
}

func (V *VPlane) Parametric(parameters ...[]float64) []float64 {
	var pars []float64
	if len(parameters) > 1 || len(parameters[0]) >= 4 {
		pars = parameters[0]

	} else {
		pars = make([]float64, 4)
	}
	pars[0] = V.Normal.At(0, 0)
	pars[1] = V.Normal.At(0, 1)
	pars[2] = V.Normal.At(0, 2)
	pars[3] = V.Normal.Dot(V.Point)
	return pars
}

//OtherAtom, given the index of one of the atoms bisected by the plane,returns the index of the other atom.
func (V *VPlane) OtherAtom(i int) int {
	if V.Atoms[0] == i {
		return V.Atoms[1]
	} else if V.Atoms[1] == i {
		return V.Atoms[0]
	}
	panic(fmt.Sprintf("Plane is not related to atom %d", i)) //I think this is fair enough. This should always be a bug in the program.
}

//DistanceInterVector  obtains the distance from a point o to the plane in the direction to the point d
//NOT in the direction of the vector d!
//It returns -1 if the vector never intersects the plane.
func (V *VPlane) DistanceInterVector(o, d *v3.Matrix) float64 {
	p := v3.Zeros(1)
	p.Sub(d, o)
	dot := V.Normal.Dot(p)
	if dot == 0 {
		return -1 //vector and plane never intersect
		//The other corner case, where the vector is contained in the plane, should not happen here!
	}
	n := V.Normal
	a := V.Point
	c := func(v *v3.Matrix, i int) float64 { return v.At(0, i) }
	t := (c(n, 0)*(c(o, 0)-c(a, 0)) + c(n, 1)*(c(o, 1)-c(a, 1)) + c(n, 2)*(c(o, 2)-c(a, 2))) / (c(n, 0)*c(d, 0) + c(n, 1)*c(d, 1) + c(n, 2)*c(d, 2))
	x := c(o, 0) + c(d, 0)*t
	y := c(o, 1) + c(d, 1)*t
	z := c(o, 2) + c(d, 2)*t
	p.Set(0, 0, x)
	p.Set(0, 1, y)
	p.Set(0, 2, z)
	p2 := v3.Zeros(1)
	p2.Sub(p, o)
	return p2.Norm(2)
}

//VPSlice contains a slice to pointers to VPlane.
type VPSlice []*VPlane

//AtomPlanes returns all the planes that bisect the atom with index i, and any other atom.
func (P VPSlice) AtomPlanes(i int) VPSlice {
	ret := make([]*VPlane, 0, len(P)/10) //the cap value is just a wild guess
	for _, v := range P {
		if v.Atoms[0] == i || v.Atoms[1] == i {
			ret = append(ret, v)
		}
	}
	return VPSlice(ret)
}

//PairPlane returns the plane bisecting the atoms with indexes i and j.
func (P VPSlice) PairPlane(i, j int) *VPlane {
	for _, v := range P {
		if (v.Atoms[0] == i && v.Atoms[1] == j) || (v.Atoms[0] == j && v.Atoms[1] == i) {
			return v
		}
	}
	return nil
}

//Determines whether 2 atoms are in contact, using the sum of their vdW radii (multiplied by an optional
//factor, 1.2 by default), as a cutoff.
func (P VPSlice) VdwContact(coords *v3.Matrix, mol chem.Atomer, i, j int, scan ...*AngleScan) bool {
	if len(scan) <= 0 {
		scan = append(scan, DefaultAngleScan()) //1.4 is the vdW radius of water

	}
	vdwsum := mol.Atom(i).Vdw + mol.Atom(j).Vdw
	scan[0].Cutoff = vdwsum*scan[0].VdwFactor + scan[0].Offset
	return P.Contact(coords, i, j, scan[0])

}

func distance(coords *v3.Matrix, i, j int) float64 {
	vi := coords.VecView(i)
	vj := coords.VecView(j)
	d := v3.Zeros(1)
	d.Sub(vi, vj)
	return d.Norm(2)
}

// Determines whether vectors i and j from coords are in contact, by checking that no plane is closer to i
// than the plane bisecting the ij vector, along some direction with an angle smaller than 90 degrees from the
// direction along the ij vector.
// This function is safe for concurrent use.
func (P VPSlice) Contact(coords *v3.Matrix, i, j int, anglest ...*AngleScan) bool {
	var angles []float64 //if given, should contain 3 values. First angle, last angle, and step (all in degrees).
	cutoff := 6.0        //some default (rather permisive) cutoff
	if len(anglest) == 0 {
		angles = []float64{1, 2} //these values ensure that only one angle 0 deg is tested.
	} else {
		angles = anglest[0].Angles
		cutoff = anglest[0].Cutoff

	}
	ref := P.PairPlane(i, j)
	if ref == nil {
		return false //the plane was never created, probably above the initial distance cutoff
	}
	if distance(coords, i, j) > anglest[0].Cutoff {
		ref.NotContact = true
		return false //profiling is good! This made everything over an order of magnitude faster xD

	}
	//no need to test things that are too far anyway.
	planes := P.AtomPlanes(i)
	notcontact := make([]bool, len(planes))
	for k, v := range planes {
		if v.Distance > anglest[0].Cutoff/2 {
			notcontact[k] = true
			//		planes[k].NotContact = true
		} else {
			notcontact[k] = false

			//		planes[k].NotContact = false
		}
	}
	ati := coords.VecView(i)
	atj := coords.VecView(j)
	for i, v := range planes {
		if v.Atoms[1] == j || v.Atoms[0] == j {
			continue //shouldn't be needed, really
		}
		if notcontact[i] {
			continue
		}
		var fname []string
		if len(anglest) > 0 && anglest[0].Test {
			fname = append(fname, fmt.Sprintf("test_%d_%d.xyz", i, j)) //test/debug
		}
		blocked := ConeBlock(ref, v, ati, atj, cutoff, angles, fname...)
		if blocked {
			ref.NotContact = true
			return false
		}
	}
	ref.Contact = true
	return true
}

//checks if the ith plane is repeated elsewhere
func (P VPSlice) IsRepeated(p *VPlane) bool {
	for _, v := range P {
		if p.Atoms[0] == v.Atoms[1] && p.Atoms[1] == v.Atoms[0] {
			return true
		}
		if p.Atoms[0] == v.Atoms[0] && p.Atoms[1] == v.Atoms[1] {
			//this should never happen!
			return true
		}
	}
	return false
}

//All contacts returns a list of pairs of atoms that are in contact
//It is assumed that the P slice doesn't contain any non-contact or any repeated plane
func (P VPSlice) AllContacts() [][2]int {
	r := make([][2]int, 0, len(P))
	for _, v := range P {
		r = append(r, [2]int{v.Atoms[0], v.Atoms[1]})
	}
	return r
}

func (P VPSlice) IsBlocked(p *VPlane, c *v3.Matrix, mustbeatom int, anglest ...*AngleScan) bool {
	var as *AngleScan
	if len(anglest) == 0 {
		as = DefaultAngleScan()
	}
	ai := c.VecView(p.Atoms[0])
	aj := c.VecView(p.Atoms[1])

	for _, v := range P {
		if v.Atoms[1] == p.Atoms[0] && v.Atoms[0] == p.Atoms[1] {
			continue //This is just the same plane
		}
		//if a "mustbeatom" is given, we only check the planes that contain that atom
		//if mustbeatom >= 0 && !(v.Atoms[1] == mustbeatom || v.Atoms[0] == mustbeatom) {
		//	continue

		if !(v.Atoms[1] == p.Atoms[0] || v.Atoms[0] == p.Atoms[0] || v.Atoms[1] == p.Atoms[1] || v.Atoms[0] == p.Atoms[1]) {
			continue //there is no point to check planes that have no atoms in common with the test one.
		}

		blocked := ConeBlock(p, v, ai, aj, as.Cutoff, as.Angles)
		if blocked {
			p.NotContact = true
			return true
		}
	}
	p.Contact = true
	return false

}

func (P VPSlice) ConfirmedContacts() VPSlice {
	ret := make([]*VPlane, 0, len(P)/10)
	for _, v := range P {
		if v.Contact {
			ret = append(ret, v)
		}

	}
	return ret
}

//returns a crappy xyz string from the slice of 1x3 vectors provided
//for testing purposes only
func writetestxyzstring(vec ...*v3.Matrix) string {
	form := fmt.Sprintf
	elem := []string{"O", "C", "N", "H", "P", "F", "I", "Br", "Zn", "Cu"}
	str := form("%d\n\n", len(vec))
	for i, v := range vec {
		str += form("%s %5.3f %5.3f %5.3f\n", elem[i], v.At(0, 0), v.At(0, 1), v.At(0, 2))
	}
	return str

}

//Test if the path between ati and the ref plane is blocked by the test plane
//it will scan cones at increasing angles from the ati-atj vector(from 0 to angles[0] degrees, where angles[0] should be <90)
//in angle[1] steps.
//This is a brute-force, very slow and clumsy system. But hey, I'm a chemist. I'll change it when a) there is a pure Go 3D-Voronoi
//library or b) I find the time to study computational geometry.
func ConeBlock(ref, test *VPlane, ati, atj *v3.Matrix, cutoff float64, angles []float64, testname ...string) bool {
	var xyz string
	refdist := ref.Distance
	axis := v3.Zeros(1)
	aux := v3.Zeros(1)
	ax2 := v3.Zeros(1)
	aux.Set(0, 1, 0) //any vector
	axis.Sub(atj, ati)
	aux.Sub(axis, aux) //just to ensure axis (defined as atj-ati) is not colinear with aux
	ax2.Cross(axis, aux)
	if angles[0] >= 90 { //this will never intersect the ref plane!
		angles[0] = 89
	}
	for ang := 0.0; ang <= angles[0]; ang += angles[1] {
		if ang == 0 {
			refdist = ref.Distance
			Dij_k := test.DistanceInterVector(ati, atj)
			if Dij_k > refdist || Dij_k < 0 {
				return false
			}
		}
		//now the cone
		dir, err := chem.RotateAbout(atj, ati, ax2, ang*chem.Deg2Rad)
		if err != nil {
			panic(err.Error()) //I am not 100% sure about panicking over this, but it seems like it has to be a bug in the caller
		}
		for rot := 0.0; rot < 360; rot += angles[1] {
			ndir, err := chem.RotateAbout(dir, ati, atj, rot)
			if err != nil {
				panic(err.Error())
			}
			if len(testname) > 0 {
				xyz += writetestxyzstring(ati, atj, ndir)
			}
			refdist = ref.DistanceInterVector(ati, ndir)
			if refdist > cutoff/2 && refdist < 0 {
				return true
			}

			Dij_k := test.DistanceInterVector(ati, ndir)
			if Dij_k > refdist || Dij_k < 0 { //a distance <0 means that the vector never intersects the plane.
				writetest(xyz, testname)
				return false
			}
		}
	}
	//for debugging
	writetest(xyz, testname)
	return true

}

//debug function, to be called on testing.
func writetest(xyz string, name []string) {
	if len(name) > 0 {
		f, err := os.Create(name[0])
		if err != nil {
			log.Println(err.Error())
		}
		f.WriteString(xyz)
		f.Close()

	}

}

//
func PlaneBetweenAtoms(at1, at2 *v3.Matrix, i, j int) *VPlane {
	ret := &VPlane{}
	ret.Atoms = []int{i, j}
	ret.Normal = v3.Zeros(1)
	ret.Normal.Sub(at2, at1)
	ret.Distance = ret.Normal.Norm(2) / 2.0
	ret.Point = v3.Zeros(1)
	ret.Point.Add(at1, at2)
	ret.Point.Scale(0.5, ret.Point) //I'm actually not 100% sure about this! D:
	return ret
}

//get all planes between all possible pairs of atoms
//which are not farther away from each other than cutoff
func GetPlanes(atoms *v3.Matrix, mol chem.Atomer, cutoff float64, noHs ...bool) VPSlice {
	var noH bool //false by default (it's zero-value)
	if len(noHs) != 0 {
		noH = noHs[0]
	}
	var ai, aj *v3.Matrix
	//ai := v3.Zeros(1)
	//aj := v3.Zeros(1)
	tot := atoms.NVecs()
	planes := make([]*VPlane, 0, combin.Binomial(tot, 2)/2) //I couldn't resist adding the binomial coefficient, sorry.
	for i := 0; i < tot; i++ {
		ai = atoms.VecView(i)
		ati := mol.Atom(i)
		if noH && ati.Symbol == "H" {
			continue
		}
		for j := i + 1; j < tot; j++ {
			aj = atoms.VecView(j)
			atj := mol.Atom(j)
			if noH && atj.Symbol == "H" {
				continue
			}
			t := PlaneBetweenAtoms(ai, aj, i, j)
			if t.Distance*2.0 > cutoff {
				continue
			}
			planes = append(planes, t)

		}
	}
	return VPSlice(planes)
}
