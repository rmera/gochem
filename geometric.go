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

import (
	"fmt"
	"math"
	"sort"

	"github.com/gonum/matrix/mat64"
	"github.com/rmera/gochem/v3"
)



//NOTE: For many of these functions we could ask for a buffer vector in the arguments in order to reduce
//memory allocation.


//Angle takes 2 vectors and calculate the angle in radians between them
//It does not check for correctness or return errors!
func Angle(v1, v2 *v3.Matrix) float64 {
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

//RotatorToNewY takes a set of coordinates (mol) and a vector (y). It returns
//a rotation matrix that, when applied to mol, will rotate it around the Z axis
//in such a way that the projection of newy in the XY plane will be aligned with
//the Y axis.
func RotatorAroundZToNewY(newy *v3.Matrix) (*v3.Matrix, error) {
	nr, nc := newy.Dims()
	if nc != 3 || nr != 1 {
		return nil, CError{"Wrong newy vector", []string{"RotatorAroundZtoNewY"}}
	}
	if nc != 3 {
		return nil, CError{"Wrong mol vector", []string{"RotatorAroundZtoNewY"}} //this one doesn't seem reachable

	}
	gamma := math.Atan2(newy.At(0, 0), newy.At(0, 1))
	singamma := math.Sin(gamma)
	cosgamma := math.Cos(gamma)
	operator := []float64{cosgamma, singamma, 0,
		-singamma, cosgamma, 0,
		0, 0, 1}
	return v3.NewMatrix(operator)

}

//RotatorAroundZ returns an operator that will rotate a set of
//coordinates by gamma radians around the z axis.
func RotatorAroundZ(gamma float64) (*v3.Matrix, error) {
	singamma := math.Sin(gamma)
	cosgamma := math.Cos(gamma)
	operator := []float64{cosgamma, singamma, 0,
		-singamma, cosgamma, 0,
		0, 0, 1}
	return v3.NewMatrix(operator)

}

//RotatorToNewZ takes a matrix a row vector (newz).
//It returns a linear operator such that, when applied to a matrix mol ( with the operator on the right side)
//it will rotate mol such that the z axis is aligned with newz.
func RotatorToNewZ(newz *v3.Matrix) *v3.Matrix {
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
	finalop, _ := v3.NewMatrix(operator) //we are hardcoding opperator so it must have the right dimensions.
	return finalop

}

//RotatorTranslatorToSuper superimposes the set of cartesian coordinates given as the rows of the matrix test on the gnOnes of the rows
//of the matrix templa. Returns the transformed matrix, the rotation matrix, 2 translation row vectors
//For the superposition plus an error. In order to perform the superposition, without using the transformed
//the first translation vector has to be added first to the moving matrix, then the rotation must be performed
//and finally the second translation has to be added.
//This is a low level function, although one can use it directly since it returns the transformed matrix.
//The math for this function is by Prof. Veronica Jimenez-Curihual, University of Concepcion, Chile.
func RotatorTranslatorToSuper(test, templa *v3.Matrix) (*v3.Matrix, *v3.Matrix, *v3.Matrix, *v3.Matrix, error) {
	tmr, tmc := templa.Dims()
	tsr, tsc := test.Dims()
	if tmr != tsr || tmc != 3 || tsc != 3 {
		return nil, nil, nil, nil, CError{"goChem: Ill-formed matrices", []string{"RotatorTranslatorToSuper"}}
	}
	var Scal float64
	Scal = float64(1.0) / float64(tmr)
	j := gnOnes(tmr, 1) //Mass is not important for this matter so we'll just use this.
	ctest, distest, err := MassCenter(test, test, j)
	if err != nil {
		return nil, nil, nil, nil, errDecorate(err, "RotatorTranslatorToSuper")
	}
	ctempla, distempla, err := MassCenter(templa, templa, j)
	if err != nil {
		return nil, nil, nil, nil, errDecorate(err, "RotatorTranslatorToSuper")

	}
	Mid := gnEye(tmr)
	jT := gnT(j)
	ScaledjProd := gnMul(j, jT)
	ScaledjProd.Scale(Scal, ScaledjProd)
	aux2 := gnMul(gnT(ctempla), Mid)
	r, _ := aux2.Dims()
	Maux := v3.Zeros(r)
	Maux.Mul(aux2, ctest)
	Maux.Tr() //Dont understand why this is needed
	factors := mat64.SVD(v3.Matrix2Dense(Maux), appzero, math.SmallestNonzeroFloat64, true, true)
	U := factors.U
	V := factors.V
	//	if err != nil {
	//		return nil, nil, nil, nil, err  //I'm not sure what err is this one
	//	}
	U.Scale(-1, U)
	V.Scale(-1, V)
	//SVD gives different results here than in numpy. U and V are multiplide by -1 in one of them
	//and gomatrix gives as V the transpose of the matrix given as V by numpy. I guess is not an
	//error, but don't know for sure.
	vtr, _ := V.Dims()
	Rotation := v3.Zeros(vtr)
	Rotation.Mul(V, gnT(U))
	Rotation.Tr() //Don't know why does this work :(
	RightHand := gnEye(3)
	if det(Rotation) < 0 {
		RightHand.Set(2, 2, -1)
		Rotation.Mul(V, RightHand)
		Rotation.Mul(Rotation, gnT(U)) //If I get this to work Ill arrange so gnT(U) is calculated once, not twice as now.
		Rotation.Tr() //TransposeTMP contains the transpose of the original Rotation      //Same, no ide why I need this
		//return nil, nil, nil, nil, fmt.Errorf("Got a reflection instead of a translations. The objects may be specular images of each others")
	}
	jT.Scale(Scal, jT)
	subtempla := v3.Zeros(tmr)
	subtempla.Copy(ctempla)
	sub := v3.Zeros(ctest.NVecs())
	sub.Mul(ctest, Rotation)
	subtempla.Sub(subtempla, sub)
	jtr, _ := jT.Dims()
	Translation := v3.Zeros(jtr)
	Translation.Mul(jT, subtempla)
	Translation.Add(Translation, distempla)
	//This alings the transformed with the original template, not the mean centrate one
	transformed := v3.Zeros(ctest.NVecs())
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
func RMSD(test, template *v3.Matrix) (float64, error) {
	//This is a VERY naive implementation.
	tmr, tmc := template.Dims()
	tsr, tsc := test.Dims()
	if tmr != tsr || tmc != 3 || tsc != 3 {
		return 0, fmt.Errorf("Ill formed matrices for RMSD calculation")
	}
	tr := tmr
	ctempla := v3.Zeros(template.NVecs())
	ctempla.Copy(template)
	//the maybe thing might not be needed since we check the dimensions before.
	f := func() { ctempla.Sub(ctempla, test) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return 0, CError{err.Error(), []string{"v3.Matrix.Sub", "RMSD"}}
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

//Dihedral calculates the dihedral between the points a, b, c, d, where the first plane
//is defined by abc and the second by bcd.
func Dihedral(a, b, c, d *v3.Matrix) float64 {
	all := []*v3.Matrix{a, b, c, d}
	for number, point := range all {
		pr, pc := point.Dims()
		if point == nil {
			panic(PanicMsg(fmt.Sprintf("goChem-Dihedral: Vector %d is nil", number)))
		}
		if pr != 1 || pc != 3 {
			panic(PanicMsg(fmt.Sprintf("goChem-Dihedral: Vector %d has invalid shape", number)))
		}
	}
	//bma=b minus a
	bma := v3.Zeros(1)
	cmb := v3.Zeros(1)
	dmc := v3.Zeros(1)
	bmascaled := v3.Zeros(1)
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
	return linear_distortion, circular_distortion, errDecorate(err, "RhoShapeIndexes")
}

//Rhos returns the semiaxis of the elipoid of inertia given the eigenvectors of the moment tensor.
func Rhos(evals []float64) ([]float64, error) {
	rhos := sort.Float64Slice{invSqrt(evals[0]), invSqrt(evals[1]), invSqrt(evals[2])}
	if evals[2] <= appzero {
		return rhos[:], CError{"goChem: Molecule colapsed to a single point. Check for blackholes", []string{"Rhos"}}
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
//Plane that best contains the molecule. Notice that the function can't possibly check
//that the vectors are sorted. The P at the end of the name is for Performance. If
//That is not an issue it is safer to use the BestPlane function that wraps this one.
func BestPlaneP(evecs *v3.Matrix) (*v3.Matrix, error) {
	evr, evc := evecs.Dims()
	if evr != 3 || evc != 3 {
		return evecs, CError{"goChem: Eigenvectors matrix must be 3x3", []string{"BestPlaneP"}} //maybe this should be a panic
	}
	v1 := evecs.VecView(2)
	v2 := evecs.VecView(1)
	normal := v3.Zeros(1)
	normal.Cross(v1, v2)
	return normal, nil
}

//BestPlane returns a row vector that is normal to the plane that best contains the molecule
//if passed a nil Masser, it will simply set all masses to 1.
func BestPlane(coords *v3.Matrix, mol Masser) (*v3.Matrix, error) {
	var err error
	var Mmass []float64
	cr, _ := coords.Dims()
	if mol != nil {
		Mmass, err = mol.Masses()
		if err != nil {
			return nil, errDecorate(err, "BestPlane")
		}
		if len(Mmass) != cr {
			return nil, CError{fmt.Sprintf("Inconsistent coordinates(%d)/atoms(%d)", len(Mmass), cr), []string{"BestPlane"}}
		}
	}
	moment, err := MomentTensor(coords, Mmass)
	if err != nil {
		return nil, errDecorate(err, "BestPlane")
	}
	evecs, _, err := v3.EigenWrap(moment, appzero)
	if err != nil {
		return nil, errDecorate(err, "BestPlane")
	}
	normal, err := BestPlaneP(evecs)
	if err != nil {
		return nil, errDecorate(err, "BestPlane")
	}
	//MomentTensor(, mass)
	return normal, err
}

//returns a flat64 slice of the size requested filed with ones
func ones(size int) []float64 {
	slice := make([]float64, size, size)
	for k, _ := range slice {
		slice[k] = 1.0
	}
	return slice
}

//CenterOfMass returns the center of mass the atoms represented by the coordinates in geometry
//and the masses in mass, and an error. If mass is nil, it calculates the geometric center
func CenterOfMass(geometry *v3.Matrix, mass *mat64.Dense) (*v3.Matrix, error) {
	if geometry == nil {
		return nil, CError{"goChem: nil matrix to get the center of mass", []string{"CenterOfMass"}}
	}
	gr, _ := geometry.Dims()
	if mass == nil { //just obtain the geometric center
		tmp := ones(gr)
		mass = mat64.NewDense(gr, 1, tmp) //gnOnes(gr, 1)
	}
	tmp2 := ones(gr)
	gnOnesvector := mat64.NewDense(1, gr, tmp2) //gnOnes(1, gr)

	ref := v3.Zeros(gr)
	ref.ScaleByCol(geometry, mass)
	ref2 := v3.Zeros(1)
	ref2.Mul(gnOnesvector, ref)
	ref2.Scale(1.0/mass.Sum(), ref2)
	return ref2, nil
}

//MassCenter centers in in the center of mass of oref. Mass must be
//A column vector. Returns the centered matrix and the displacement matrix.
func MassCenter(in, oref *v3.Matrix, mass *mat64.Dense) (*v3.Matrix, *v3.Matrix, error) {
	or, _ := oref.Dims()
	ir, _ := in.Dims()
	if mass == nil { //just obtain the geometric center
		tmp := ones(or)
		mass = mat64.NewDense(or, 1, tmp) //gnOnes(or, 1)
	}
	ref := v3.Zeros(or)
	ref.Copy(oref)
	gnOnesvector := gnOnes(1, or)
	f := func() { ref.ScaleByCol(ref, mass) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return nil, nil, CError{err.Error(), []string{"v3.Matrix.ScaleByCol", "MassCenter"}}
	}
	ref2 := v3.Zeros(1)
	g := func() { ref2.Mul(gnOnesvector, ref) }
	if err := gnMaybe(gnPanicker(g)); err != nil {
		return nil, nil, CError{err.Error(), []string{"v3.gOnesVector", "MassCenter"}}
	}
	ref2.Scale(1.0/mass.Sum(), ref2)
	returned := v3.Zeros(ir)
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
func MomentTensor(A *v3.Matrix, massslice []float64) (*v3.Matrix, error) {
	ar, ac := A.Dims()
	var err error
	var mass *mat64.Dense
	if massslice == nil {
		mass = gnOnes(ar, 1)
	} else {
		mass = mat64.NewDense(ar, 1, massslice)
		//		if err != nil {
		//			return nil, err
		//		}
	}
	center, _, err := MassCenter(A, v3.Dense2Matrix(gnCopy(A)), mass)
	if err != nil {
		return nil, errDecorate(err, "MomentTensor")
	}
	sqrmass := gnZeros(ar, 1)
	//	sqrmass.Pow(mass,0.5)
	pow(mass, sqrmass, 0.5) //the result is stored in sqrmass
	//	fmt.Println(center,sqrmass) ////////////////////////
	center.ScaleByCol(center, sqrmass)
	//	fmt.Println(center,"scaled center")
	centerT := gnZeros(ac, ar)
	centerT.Copy(center.T())
	moment := gnMul(centerT, center)
	return v3.Dense2Matrix(moment), nil
}

//Projection returns the projection of test in ref.
func Projection(test, ref *v3.Matrix) *v3.Matrix {
	rr, _ := ref.Dims()
	Uref := v3.Zeros(rr)
	Uref.Unit(ref)
	scalar := test.Dot(Uref) //math.Abs(la)*math.Cos(angle)
	Uref.Scale(scalar, Uref)
	return Uref
}

//AntiProjection returns a vector in the direction of ref with the magnitude of
//a vector A would have if |test| was the magnitude of its projection
//in the direction of test.
func AntiProjection(test, ref *v3.Matrix) *v3.Matrix{
	rr,_:=ref.Dims()
	testnorm:=test.Norm(0)
	Uref:=v3.Zeros(rr)
	Uref.Unit(ref)
	scalar:=test.Dot(Uref)
	scalar=(testnorm*testnorm)/scalar
	Uref.Scale(scalar,Uref)
	return Uref
}

//SelCone, Given a set of cartesian points in sellist, obtains a vector "plane" normal to the best plane passing through the points.
//It selects atoms from the set A that are inside a cone in the direction of "plane" that starts from the geometric center of the cartesian points,
//and has an angle of angle (radians), up to a distance distance. The cone is approximated by a set of radius-increasing cilinders with height thickness.
//If one starts from one given point, 2 cgnOnes, one in each direction, are possible. If whatcone is 0, both cgnOnes are considered.
//if whatcone<0, only the cone opposite to the plane vector direction. If whatcone>0, only the cone in the plane vector direction.
//the 'initial' argument  allows the construction of a truncate cone with a radius of initial.
func SelCone(B, selection *v3.Matrix, angle, distance, thickness, initial float64, whatcone int) []int {
	A := v3.Zeros(B.NVecs())
	A.Copy(B) //We will be altering the input so its better to work with a copy.
	ar, _ := A.Dims()
	selected := make([]int, 0, 3)
	neverselected := make([]int, 0, 30000)     //waters that are too far to ever be selected
	nevercutoff := distance / math.Cos(angle)  //cutoff to be added to neverselected
	A, _, err := MassCenter(A, selection, nil) //Centrate A in the geometric center of the selection, Its easier for the following calculations
	if err != nil {
		panic(PanicMsg(err.Error()))
	}
	selection, _, _ = MassCenter(selection, selection, nil) //Centrate the selection as well
	plane, err := BestPlane(selection, nil)                 //I have NO idea which direction will this vector point. We might need its negative.
	if err != nil {
		panic(PanicMsg(err.Error()))
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
