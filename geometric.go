/*
 * geometric.go, part of gochem
 *
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 *
*/
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

package chem

import "fmt"
import "math"
import "sort"

/*
   """projects any vector V in a plane with contains the 0 and which normal vector is  N. V and N must be numeric arrays"""
   def __init__(self,N): #set the operator
           self.file1=[1-(pow(N[0],2)/pow(norm(N),2)),(N[0]*N[1])/pow(norm(N),2),(N[0]*N[2])/pow(norm(N),2)]
           self.file2=[(N[1]*N[0])/pow(norm(N),2),1-(pow(N[1],2)/pow(norm(N),2)),(N[1]*N[2])/pow(norm(N),2)]
           self.file3=[(N[2]*N[0])/pow(norm(N),2),(N[2]*N[1])/pow(norm(N),2),1-(pow(N[2],2)/pow(norm(N),2))]
           self.op=array([self.file1,self.file2,self.file3])
   def project(self,V):
           self.projection=dot(self.op,V)
           return self.projection


*/
//func GetShadow()

//Angle takes 2 vectors and calculate the angle in radians between them
//It does not check for correctness or return errors!
func Angle(v1, v2 *VecMatrix) float64 {
	normproduct := v1.Norm(0) * v2.Norm(0)
	dotprod := v1.Dot(v2)
	argument := dotprod / normproduct
	//Take care of floating point math errors
	if math.Abs(argument-1) <= appzero {
		argument = 1
	} else if math.Abs(argument+1) <= appzero {
		argument = -1
	}
	//fmt.Println(dotprod/normproduct,argument) //dotprod/normproduct, dotprod, normproduct,v1.TwoNorm(),v2.TwoNorm())
	angle := math.Acos(argument)
	if math.Abs(angle) <= appzero {
		return 0.00
	}
	return angle
}

//GetRotateToNewY takes a set of coordinates (mol) and a vector (y). It returns
//a rotation matrix that, when applied to mol, will rotate it around the Z axis
//in such a way that the projection of newy in the XY plane will be aligned with
//the Y axis.
func RotatorAroundZToNewY(newy *VecMatrix) (*VecMatrix, error) {
	nr, nc := newy.Dims()
	if nc != 3 || nr != 1 {
		return nil, fmt.Errorf("Wrong newy vector")
	}
	if nc != 3 {
		return nil, fmt.Errorf("Wrong mol vector")
	}
	gamma := math.Atan2(newy.At(0, 0), newy.At(0, 1))
	singamma := math.Sin(gamma)
	cosgamma := math.Cos(gamma)
	operator := []float64{cosgamma, singamma, 0,
		-singamma, cosgamma, 0,
		0, 0, 1}
	return NewVecs(operator)

}

//GetRotateAroundZ returns an operator that will rotate a set of
//coordinates by gamma radians around the z axis.
func RotatorAroundZ(gamma float64) (*VecMatrix, error) {
	singamma := math.Sin(gamma)
	cosgamma := math.Cos(gamma)
	operator := []float64{cosgamma, singamma, 0,
		-singamma, cosgamma, 0,
		0, 0, 1}
	return NewVecs(operator)

}

//RotatorToNewZ takes a matrix a row vector (newz).
//It returns a linear operator such that, when applied to a matrix mol ( with the operator on the right side)
//it will rotate mol such that the z axis is aligned with newz.
func RotatorToNewZ(newz *VecMatrix) *VecMatrix {
	r, c := newz.Dims()
	if c != 3 || r != 1 {
		panic("Wrong newz vector")
	}
	normxy := math.Sqrt(math.Pow(newz.At(0, 0), 2) + math.Pow(newz.At(0, 1), 2))
	theta := math.Atan2(normxy, newz.At(0, 2))      //Around the new y
	phi := math.Atan2(newz.At(0, 1), newz.At(0, 0)) //First around z
	psi := 0.000000000000                           // second around z
	sinphi := math.Sin(phi)
	cosphi := math.Cos(phi)
	sintheta := math.Sin(theta)
	costheta := math.Cos(theta)
	sinpsi := math.Sin(psi)
	cospsi := math.Cos(psi)
	operator := []float64{cosphi*costheta*cospsi - sinphi*sinpsi, -sinphi*cospsi - cosphi*costheta*sinpsi, cosphi * sintheta,
		sinphi*costheta*cospsi + cosphi*sinpsi, -sinphi*costheta*sinpsi + cosphi*cospsi, sintheta * sinphi,
		-sintheta * cospsi, sintheta * sinpsi, costheta}
	finalop, _ := NewVecs(operator) //we are hardcoding opperator so it must have the right dimensions.
	return finalop

}

//RotatorTranslatorToSuper superimposes the set of cartesian coordinates given as the rows of the matrix test on the gnOnes of the rows
//of the matrix templa. Returns the transformed matrix, the rotation matrix, 2 translation row vectors
//For the superposition plus an error. In order to perform the superposition, without using the transformed
//the first translation vector has to be added first to the moving matrix, then the rotation must be performed
//and finally the second translation has to be added.
//This is a low level function, although one can use it directly since it returns the transformed matrix.
//The math for this function is by Prof. Veronica Jimenez-Curihual, University of Concepcion, Chile.
func RotatorTranslatorToSuper(test, templa *VecMatrix) (*VecMatrix, *VecMatrix, *VecMatrix, *VecMatrix, error) {
	tmr, tmc := templa.Dims()
	tsr, tsc := test.Dims()
	if tmr != tsr || tmc != 3 || tsc != 3 {
		return nil, nil, nil, nil, fmt.Errorf("GetSuper: Ill-formed matrices")
	}
	var Scal float64
	Scal = float64(1.0) / float64(tmr)
	j := gnOnes(tmr, 1) //Mass is not important for this matter so we'll just use this.
	ctest, distest, err := MassCentrate(test, test, j)
	if err != nil {
		return nil, nil, nil, nil, err
	}
	ctempla, distempla, err := MassCentrate(templa, templa, j)
	if err != nil {
		return nil, nil, nil, nil, err
	}
	Mid := gnEye(tmr)
	jT := gnT(j)
	ScaledjProd := gnMul(j, jT)
	ScaledjProd.Scale(Scal, ScaledjProd)
	aux2 := gnMul(gnT(ctempla), Mid)
	r, _ := aux2.Dims()
	Maux := ZeroVecs(r)
	Maux.Mul(aux2, ctest)
	Maux.TCopy(Maux) //Dont understand why this is needed
	U, _, Vt := gnSVD(vecMatrix2chemDense(Maux))
	if err != nil {
		return nil, nil, nil, nil, err
	}
	U.Scale(-1, U)
	Vt.Scale(-1, Vt)
	//SVD gives different results here than in numpy. U and Vt are multiplide by -1 in one of them
	//and gomatrix gives as Vt the transpose of the matrix given as Vt by numpy. I guess is not an
	//error, but don't know for sure.
	vtr, _ := Vt.Dims()
	Rotation := ZeroVecs(vtr)
	Rotation.Mul(Vt, gnT(U))
	Rotation.TCopy(Rotation) //Don't know why does this work :(
	if det(Rotation) < 0 {
		return nil, nil, nil, nil, fmt.Errorf("Got a reflection instead of a translations. The objects may be specular images of each others")
	}
	jT.Scale(Scal, jT)
	subtempla := ZeroVecs(tmr)
	subtempla.Copy(ctempla)
	sub := ZeroVecs(ctest.NVecs())
	sub.Mul(ctest, Rotation)
	subtempla.Sub(subtempla, sub)
	jtr, _ := jT.Dims()
	Translation := ZeroVecs(jtr)
	Translation.Mul(jT, subtempla)
	Translation.Add(Translation, distempla)
	//This allings the transformed with the original template, not the mean centrate one
	transformed := ZeroVecs(ctest.NVecs())
	transformed.Mul(ctest, Rotation)
	transformed.AddVec(transformed, Translation)
	//end transformed
	distest.Scale(-1, distest)
	return transformed, Rotation, distest, Translation, nil
}

/*
//I keep this just in case I manage to fix it at some point
func rmsd_fail(test, template *matrix.DenseMatrix) (float64, error) {
	ctempla := template.Copy()
	err := ctempla.Subtract(test)
	if err != nil {
		return 9999999, err
	}
	dev := matrix.ParallelProduct(ctempla.Transpose(), ctempla)
	RMSDv := dev.Trace()
	RMSDv = math.Sqrt(RMSDv)
	return RMSDv, nil
}
*/

//RMSD returns the RSMD (root of the mean square deviation) for the sets of cartesian
//coordinates in test and template.
func RMSD(test, template *VecMatrix) (float64, error) {
	//This is a VERY naive implementation.
	tmr, tmc := template.Dims()
	tsr, tsc := test.Dims()
	if tmr != tsr || tmc != 3 || tsc != 3 {
		return 0, fmt.Errorf("Ill formed matrices for RMSD calculation")
	}
	tr := tmr
	ctempla := ZeroVecs(template.NVecs())
	ctempla.Copy(template)
	//the maybe thing might not be needed since we check the dimensions before.
	f := func() { ctempla.Sub(ctempla, test) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return 0, err
	}
	var RMSD float64
	for i := 0; i < template.NVecs(); i++ {
		temp := ctempla.VecView(i)
		RMSD += math.Pow(temp.Norm(0), 2)
	}
	RMSD = RMSD / float64(tr)
	RMSD = math.Sqrt(RMSD)
	return RMSD, nil
}

//Dihedral calculate the dihedral between the points a, b, c, d, where the first plane
//is defined by abc and the second by bcd.
func Dihedral(a, b, c, d *VecMatrix) float64 {
	all := []*VecMatrix{a, b, c, d}
	for number, point := range all {
		pr, pc := point.Dims()
		if point == nil {
			panic(fmt.Sprintf("Vector %d is nil", number))
		}
		if pr != 1 || pc != 3 {
			panic(fmt.Sprintf("Vector %d has invalid shape", number))
		}
	}
	//bma=b minus a
	bma := ZeroVecs(1)
	cmb := ZeroVecs(1)
	dmc := ZeroVecs(1)
	bmascaled := ZeroVecs(1)
	bma.Sub(b, a)
	cmb.Sub(c, b)
	dmc.Sub(d, c)
	bmascaled.Scale(cmb.Norm(0), bma)
	first := bmascaled.Dot(cross(cmb, dmc))
	v1 := cross(bma, cmb)
	v2 := cross(cmb, dmc)
	second := v1.Dot(v2)
	dihedral := math.Atan2(first, second)
	return dihedral
}

/***Shape indicator functions***/
//const appzero float64 = 0.0000001 //used to correct floating point
//errors. Everything equal or less than this is considered zero.

//point comparisons

//RhoShapeIndexes Get shape indices based on the axes of the elipsoid of inertia.
//Based on the work of Taylor et al., .(1983), J Mol Graph, 1, 30
//This function has NOT been tested.
func RhoShapeIndexes(evals []float64) (float64, float64, error) {
	rhos, err := Rhos(evals)
	linear_distortion := (1 - (rhos[1] / rhos[0])) * 100   //Prolate
	circular_distortion := (1 - (rhos[2] / rhos[0])) * 100 //Oblate
	return linear_distortion, circular_distortion, err
}

//Rhos returns the semiaxis of the elipoid of inertia given the eigenvectors of the moment tensor.
func Rhos(evals []float64) ([]float64, error) {
	rhos := sort.Float64Slice{invSqrt(evals[0]), invSqrt(evals[1]), invSqrt(evals[2])}
	if evals[2] <= appzero {
		return rhos[:], fmt.Errorf("Molecule colapsed to a single point. Check for blackholes.")
	}
	sort.Sort(rhos[:])
	//This loop reversing loop is almost verbatin from Effective Go
	//(http://golang.org/doc/effective_go.html)
	for i, j := 0, len(rhos)-1; i < j; i, j = i+1, j-1 {
		rhos[i], rhos[j] = rhos[j], rhos[i]
	}
	return rhos[:], nil
}

/**Other geometrical**/

//BestPlaneP takes sorted evecs, according to the eval,s and returns a row vector that is normal to the
//Plane that best contains the molecule. Note that the function can't possibly check
//That the vectors are sorted!. The P at the end of the name is for Performance. If
//That is not an issue it is safer to use the BestPlane function that wraps this one.
func BestPlaneP(evecs *VecMatrix) (*VecMatrix, error) {
	evr, evc := evecs.Dims()
	if evr != 3 || evc != 3 {
		return evecs, fmt.Errorf("Eigenvectors matrix must be 3x3")
	}
	v1 := evecs.VecView(2)
	v2 := evecs.VecView(1)
	normal := ZeroVecs(1)
	normal.Cross(v1, v2)
	return normal, nil
}

//BestPlane returns a row vector that is normal to the plane that best contains the molecule
//if passed a nil Ref, it will simply set all masses to 1.
func BestPlane(coords *VecMatrix, mol ReadRef) (*VecMatrix, error) {
	var err error
	var Mmass []float64
	cr, _ := coords.Dims()
	if mol != nil {
		if mol.Len() != cr {
			return nil, fmt.Errorf("Inconsistent coordinates(%d)/atoms(%d)", mol.Len(), cr)
		}
		Mmass, err = mol.Masses()
		if err != nil {
			return nil, err
		}
	} else {
		//Mmass=matrix.gnOnes(coords.Rows(),1)
	}
	moment, err := MomentTensor(coords, Mmass)
	if err != nil {
		return nil, err
	}
	evecs, _, err := gnEigen(moment, appzero)
	if err != nil {
		return nil, err
	}
	normal, err := BestPlaneP(evecs)
	if err != nil {
		return nil, err
	}
	//MomentTensor(, mass)
	return normal, err
}

//CenterOfMass returns the center of mass the atoms represented by the coordinates in geometry
//and the masses in mass, and an error. If mass is nil, it calculates the geometric center
func CenterOfMass(geometry *VecMatrix, mass *chemDense) (*VecMatrix, error) {
	if geometry == nil {
		return nil, fmt.Errorf("nil matrix to get the center of mass")
	}
	gr, _ := geometry.Dims()
	if mass == nil { //just obtain the geometric center
		mass = gnOnes(gr, 1)
	}
	gnOnesvector := gnOnes(1, gr)
	ref := ZeroVecs(gr)
	ref.ScaleByCol(geometry, mass)
	ref2 := ZeroVecs(1)
	ref2.Mul(gnOnesvector, ref)
	ref2.Scale(1.0/mass.Sum(), ref2)
	return ref2, nil
}

//MassCentrate centers in in the center of mass of oref. Mass must be
//A column vector. Returns the centered matrix and the displacement matrix.
func MassCentrate(in, oref *VecMatrix, mass *chemDense) (*VecMatrix, *VecMatrix, error) {
	or, _ := oref.Dims()
	ir, _ := in.Dims()
	if mass == nil { //just obtain the geometric center
		mass = gnOnes(or, 1)
	}
	ref := ZeroVecs(or)
	ref.Copy(oref)
	gnOnesvector := gnOnes(1, or)
	f := func() { ref.ScaleByCol(ref, mass) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, nil, err
	}
	ref2 := ZeroVecs(1)
	g := func() { ref2.Mul(gnOnesvector, ref) }
	if err := gnMaybe(gnPanicker(g)); err != nil {
		return nil, nil, err
	}
	ref2.Scale(1.0/mass.Sum(), ref2)
	returned := ZeroVecs(ir)
	returned.Copy(in)
	returned.SubVec(returned, ref2)
	/*	for i := 0; i < ir; i++ {
			if err := returned.GetRowVector(i).Subtract(ref2); err != nil {
				return nil, nil, err
			}
		}
	*/
	return returned, ref2, nil
}

//MomentTensor returns the moment tensor for a matrix A of coordinates and a column
//vector mass with the respective massess.
func MomentTensor(A *VecMatrix, massslice []float64) (*VecMatrix, error) {
	ar, ac := A.Dims()
	var err error
	var mass *chemDense
	if massslice == nil {
		mass = gnOnes(ar, 1)
	} else {
		mass, err = newchemDense(massslice, ar, 1)
		if err != nil {
			return nil, err
		}
	}
	center, _, err := MassCentrate(A, chemDense2VecMatrix(gnCopy(A)), mass)
	if err != nil {
		return nil, err
	}
	sqrmass := gnZeros(ar, 1)
	sqrmass.Pow(mass, 0.5)
	//	fmt.Println(center,sqrmass) ////////////////////////
	center.ScaleByCol(center, sqrmass)
	//	fmt.Println(center,"scaled center")
	centerT := gnZeros(ac, ar)
	centerT.TCopy(center)
	moment := gnMul(centerT, center)
	return chemDense2VecMatrix(moment), err
}

//The projection of test in ref.
func Projection(test, ref *VecMatrix) *VecMatrix {
	rr, _ := ref.Dims()
	Uref := ZeroVecs(rr)
	Uref.Unit(ref)
	scalar := test.Dot(Uref) //math.Abs(la)*math.Cos(angle)
	Uref.Scale(scalar, Uref)
	return Uref
}

//Given a set of cartesian points in sellist, obtains a vector "plane" normal to the best plane passing through the points.
//It selects atoms from the set A that are inside a cone in the direction of "plane" that starts from the geometric center of the cartesian points,
//and has an angle of angle (radians), up to a distance distance. The cone is approximated by a set of radius-increasing cilinders with height thickness.
//If one starts from one given point, 2 cgnOnes, one in each direction, are possible. If whatcone is 0, both cgnOnes are considered.
//if whatcone<0, only the cone opposite to the plane vector direction. If whatcone>0, only the cone in the plane vector direction.
//the 'initial' argument  allows the construction of a truncate cone with a radius of initial.
func SelCone(B, selection *VecMatrix, angle, distance, thickness, initial float64, whatcone int) []int {
	A := ZeroVecs(B.NVecs())
	A.Copy(B) //We will be altering the input so its better to work with a copy.
	ar, _ := A.Dims()
	selected := make([]int, 0, 3)
	neverselected := make([]int, 0, 30000)       //waters that are too far to ever be selected
	nevercutoff := distance / math.Cos(angle)    //cutoff to be added to neverselected
	A, _, err := MassCentrate(A, selection, nil) //Centrate A in the geometric center of the selection, Its easier for the following calculations
	if err != nil {
		panic(err.Error())
	}
	selection, _, _ = MassCentrate(selection, selection, nil) //Centrate the selection as well
	plane, err := BestPlane(selection, nil)                   //I have NO idea which direction will this vector point. We might need its negative.
	if err != nil {
		panic(err.Error())
	}
	for i := thickness / 2; i <= distance; i += thickness {
		maxdist := math.Tan(angle)*i + initial //this should give me the radius of the cone at this point
		for j := 0; j < ar; j++ {
			if isInInt(selected, j) || isInInt(neverselected, j) { //we dont scan things that we have already selected, or are too far
				continue
			}
			atom := A.VecView(j)
			proj := Projection(atom, plane)
			norm := proj.Norm(0)
			//Now at what side of the plane is the atom?
			angle := Angle(atom, plane)
			if whatcone > 0 {
				if angle > math.Pi/2 {
					continue
				}
			} else if whatcone < 0 {
				if angle < math.Pi/2 {
					continue
				}
			}
			if norm > i+(thickness/2.0) || norm < (i-thickness/2.0) {
				continue
			}
			proj.Sub(proj, atom)
			projnorm := proj.Norm(0)
			if projnorm <= maxdist {
				selected = append(selected, j)
			}
			if projnorm >= nevercutoff {
				neverselected = append(neverselected, j)
			}
		}
	}
	return selected
}
