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

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

const ( ///the nex
	defWaterOffset float64 = 2.8 //1 water radius added to the radii of each atom in contact, so 2 water radii.
	defOffset      float64 = 0.0
	defVdwFactor   float64 = 1.2
	defMaxAngle    float64 = 87
	defAngleStep   float64 = 5
	defCutoff      float64 = 5 //rather permissive
)

//ScanOptions contains options to perform angle scans to see if there is an angle in which 2 atoms
//are in direct contact (i.e. if part of the plane bisecting them is part of the Voronoi polihedra for the system).
type ScanOptions struct {
	NoH           bool
	Subset        []int
	VdwFactor     float64
	vdwSum        float64
	Offset        float64
	WaterContacts bool
	Angles        []float64 //last angle, step between angles, in degrees. A full scan would be 0 to 90
	Cutoff        float64   //distance cutoff, if the vectors at the given angle are farther than this, the angle is ignored.
}

//DefaultScanOptions returns the default setting for an ScanOptions
func DefaultScanOptions() *ScanOptions {
	return &ScanOptions{Offset: defOffset, VdwFactor: defVdwFactor, Angles: []float64{defMaxAngle, defAngleStep}, Cutoff: defCutoff}
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

//DistanceInterVector  obtains the distance from a point o to the plane in the direction to the point d
//NOT in the direction of the vector d!
//It returns -1 if the vector never intersects the plane.
//I think now this should be an average betweent eh posiotion of o and d.
func (V *VPlane) DistanceInterVector(o, d *v3.Matrix) float64 {
	od := v3.Zeros(1)
	od.Sub(d, o)
	od.Unit(od)
	dot := V.Normal.Dot(od)
	if dot == 0 {
		return -1 //odtor and plane never intersect
		//The other corner case, where the odtor is contained in the plane, should not happen,
		//but it would also lead to a dot==0
	}
	P := V.Parametric()
	subterm := P[0]*o.At(0, 0) + P[1]*o.At(0, 1) + P[2]*o.At(0, 2)
	deno := P[0]*od.At(0, 0) + P[1]*od.At(0, 1) + P[2]*od.At(0, 2)
	m := (P[3] - subterm) / deno
	//	if m < 0 {
	//		log.Println("m<0!", m)
	//	}
	return m //The sign of this number is pretty important. if your have this situation:  / .   / If you don't consider the sign
	//it could seem that the leftmost plane is blocking the way between the point and the rightmost one, as it is closer,
	//but in the wrong direction.
}

//obtains the plane equation Ax+By+Cz = D for our plane, returns A, B, C and D
func (V *VPlane) Parametric(parameters ...[]float64) []float64 {
	var pars []float64
	if len(parameters) >= 1 && len(parameters[0]) >= 4 {
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
		if v.Atoms[0] == i && v.Atoms[1] == j {
			return v
		}
		//We always return the plan in such a way that it's "from the point of view" of the first atom
		if v.Atoms[0] == j && v.Atoms[1] == i {
			v.Atoms[0], v.Atoms[1] = v.Atoms[1], v.Atoms[0]
			v.Normal.Scale(-1, v.Normal)
			return v
		}

	}
	return nil
}

func distance(coords *v3.Matrix, i, j int) float64 {
	vi := coords.VecView(i)
	vj := coords.VecView(j)
	d := v3.Zeros(1)
	d.Sub(vi, vj)
	return d.Norm(2)
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

func (P VPSlice) IsBlocked(p *VPlane, c *v3.Matrix, mustbeatom int, anglest ...*ScanOptions) bool {
	var as *ScanOptions
	if len(anglest) == 0 {
		as = DefaultScanOptions()
	}
	ai := c.VecView(p.Atoms[0])
	aj := c.VecView(p.Atoms[1])
	blocked := ConeBlockSlice(p, P, ai, aj, as.Cutoff, as.Angles)
	if blocked {
		p.NotContact = true
		return true
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

//Test if the path between ati and the ref plane is blocked by the test plane
//it will scan cones at increasing angles from the ati-atj vector(from 0 to angles[0] degrees, where angles[0] should be <90)
//in angle[1] steps.
//This is a brute-force, very slow and clumsy system. But hey, I'm a chemist. I'll change it when a) there is a pure Go 3D-Voronoi
//library or b) I find the time to study computational geometry.
func ConeBlockSlice(ref *VPlane, tests VPSlice, ati, atj *v3.Matrix, cutoff float64, angles []float64, testname ...string) bool {
	refdist := ref.Distance
	//	if refdist != ref.DistanceInterVector(ati, atj) {
	//		println("ctm", refdist, ref.DistanceInterVector(ati, atj)) ///////////
	//	}
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
	truetests := make([]*VPlane, 0, 100)
	//	fmt.Println(ref.Atoms) ////////
	//we first go through the planes to discard the ones we don't consider here
	// we collect a subset of the planes that will actually be tested
	for _, test := range tests {
		if ref.Atoms[1] == test.Atoms[0] && ref.Atoms[0] == test.Atoms[1] {
			continue //This is just the same plane
		}
		if ref.Atoms[1] == test.Atoms[1] && ref.Atoms[0] == test.Atoms[1] {
			//	println("this shouldn't happen!") /////////////
			continue //This is just the same plane
		}

		if !(ref.Atoms[1] == test.Atoms[0] || ref.Atoms[1] == test.Atoms[1] || ref.Atoms[0] == test.Atoms[0] || ref.Atoms[0] == test.Atoms[1]) {
			continue //there is no point to check planes that have no atoms in common with the test one.
		}
		truetests = append(truetests, test)

	}
	for ang := 0.0; ang <= angles[0]; ang += angles[1] {
		//now the cone
		dir, err := chem.RotateAbout(atj, ati, ax2, ang*chem.Deg2Rad)
		if err != nil {
			panic(err.Error()) //I am not 100% sure about panicking over this, but it seems like it has to be a bug in the caller
		}
		for rot := 0.0; rot < 360; rot += angles[1] {
			ndir, err := chem.RotateAbout(dir, ati, atj, rot*chem.Deg2Rad)
			if err != nil {
				panic(err.Error())
			}
			refdist = ref.DistanceInterVector(ati, ndir)
			//means at this angle there is never an intersection
			if refdist < 0 || refdist > cutoff {
				continue
			}
			blocked := false
			for _, test := range truetests {
				Dij_k := test.DistanceInterVector(ati, ndir)
				if Dij_k <= refdist && Dij_k >= 0 { //a D_ij_K <0 ==>  vector never intersects the plane.
					//				if ang == 0 && rot == 0 {
					//					fmt.Println(ref.Atoms, test.Atoms, refdist, Dij_k, ang, rot) /////////
					//				}
					blocked = true
					break
				}
			}
			if !blocked {
				return false
			}
		}
	}
	return true

}

func PlaneBetweenAtoms(at1, at2 *v3.Matrix, i, j int) *VPlane {
	ret := &VPlane{}
	ret.Atoms = []int{i, j}
	ret.Normal = v3.Zeros(1)
	ret.Normal.Sub(at1, at2)
	ret.Distance = ret.Normal.Norm(2) / 2.0
	ret.Point = v3.Zeros(1)
	ret.Point.Add(at1, at2)
	ret.Point.Scale(0.5, ret.Point) //It's the geometric center of the 2 points, so it should definitely be part of the plane.
	//println(ret.Distance, ret.DistanceInterVector(at1, at2), "SHOULD BE THE SAME!!!") ///////////
	return ret
}
