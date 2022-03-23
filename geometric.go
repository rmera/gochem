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

package chem

import (
	"fmt"
	"math"
	"sort"

	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

//NOTE: For many of these functions we could ask for a buffer vector in the arguments in order to reduce
//memory allocation.

//Angle takes 2 vectors and calculate the angle in radians between them
//It does not check for correctness or return errors!
func Angle(v1, v2 *v3.Matrix) float64 {
	normproduct := v1.Norm(2) * v2.Norm(2)
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

//RotatorAroundZToNewY takes a set of coordinates (mol) and a vector (y). It returns
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
//The math for this function is by Prof. Veronica Jimenez-Curihual, UNAB, Chile.
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
	Maux = Maux.TrRet() //Dont understand why this is needed
	svd := new(mat.SVD)
	if ok := svd.Factorize(v3.Matrix2Dense(Maux), mat.SVDFull); !ok {
		return nil, nil, nil, nil, errDecorate(fmt.Errorf("mat.SVD failed"), "RotatorTranslatorToSuper")
	}

	//matrix Maux dimensions must be 3x3
	Ud := mat.NewDense(3, 3, make([]float64, 9))
	Vd := mat.NewDense(3, 3, make([]float64, 9))
	svd.UTo(Ud)
	svd.VTo(Vd)
	U := v3.Dense2Matrix(Ud) //These will panic if the matrices were not 3x3
	V := v3.Dense2Matrix(Vd)
	/***** OLD
	factors := mat.SVD(v3.Matrix2Dense(Maux), appzero, math.SmallestNonzeroFloat64, true, true)
	U := factors.U
	V := factors.V
	****/
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
	Rotation = Rotation.TrRet() //Don't know why does this work :(
	RightHand := gnEye(3)
	if det(Rotation) < 0 {
		RightHand.Set(2, 2, -1)
		Rotation.Mul(V, RightHand)
		Rotation.Mul(Rotation, gnT(U)) //If I get this to work Ill arrange so gnT(U) is calculated once, not twice as now.
		Rotation.Tr()                  //TransposeTMP contains the transpose of the original Rotation      //Same, no ide why I need this
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

//RMSD calculates the RMSD between test and template, considering only the atoms
//present in the slices of int slices indexes. The first indexes slices will
//be assumed to contain test indexes and the second, template indexes.
//If you give only one, it will be assumed to correspondo to whatever molecule
//that has more atoms than the elements in the slice. The same number of atoms
//has to be considered for superposition in both systems.
//The objects are not superimposed before the calculation.
func RMSD(test, templa *v3.Matrix, indexes ...[]int) (float64, error) {
	var L int
	if len(indexes) == 0 || indexes[0] == nil || len(indexes[0]) == 0 {
		L = test.NVecs()
	} else {
		L = len(indexes[0])
	}
	tmp := v3.Zeros(L)
	//We don't test anything in-house. All the testing in done by MemRMSD
	rmsd, err := MemRMSD(test, templa, tmp, indexes...)
	return rmsd, err
}

//MemRMSD calculates the RMSD between test and template, considering only the atoms
//present in the slices of int slices indexes. The first indexes slices will
//be assumed to contain test indexes and the second, template indexes.
//If you give only one (it must be the first one), it will be assumed to correspond to whatever molecule
//that has more atoms than the elements in the slice.  Giving a nil or 0-lenght first slice and a non-nil second
//slice will cause MemRMSD to not consider neither of them.
//The same number of atoms
//has to be considered for the calculation in both systems.
//It does not superimpose the objects.
//To save memory, it asks for the temporary matrix it needs to be supplied:
//tmp must be Nx3 where N is the number
//of elements in testlst and templalst
func MemRMSD(test, templa, tmp *v3.Matrix, indexes ...[]int) (float64, error) {
	var ctest *v3.Matrix
	var ctempla *v3.Matrix
	if len(indexes) == 0 || indexes[0] == nil || len(indexes[0]) == 0 {
		ctest = test
		ctempla = templa
	} else if len(indexes) == 1 {
		if test.NVecs() > len(indexes[0]) {
			ctest = v3.Zeros(len(indexes[0]))
			ctest.SomeVecs(test, indexes[0])
			ctempla = templa
		} else if templa.NVecs() > len(indexes[0]) {
			ctempla = v3.Zeros(len(indexes[0]))
			ctempla.SomeVecs(templa, indexes[0])
		} else {
			return -1, fmt.Errorf("chem.memRMSD: Indexes don't match molecules")
		}
	} else {
		ctest = v3.Zeros(len(indexes[0]))
		ctest.SomeVecs(test, indexes[0])
		ctempla = v3.Zeros(len(indexes[1]))
		ctempla.SomeVecs(templa, indexes[1])
	}
	if ctest.NVecs() != ctempla.NVecs() || tmp.NVecs() != ctest.NVecs() {
		return -1, fmt.Errorf("chem.memRMSD: Ill formed matrices for memRMSD calculation")
	}
	tmp.Sub(ctest, ctempla)
	rmsd := tmp.Norm(2)
	return rmsd / math.Sqrt(float64(ctest.NVecs())), nil

}

//MemMSD calculates the MSDs between test and templa, considering only the atoms
//present in the slices of int slices indexes. If given. If only one set of indexes
//is given, it will be assumed to beling to test, if it has less elements than that
//system, or to templa,otherwise.
func MemPerAtomRMSD(test, templa, ctest, ctempla, tmp *v3.Matrix, indexes ...[]int) ([]float64, error) {
	if len(indexes) == 0 || indexes[0] == nil || len(indexes[0]) == 0 {
		ctest = test
		ctempla = templa
	} else if len(indexes) == 1 {
		if test.NVecs() > len(indexes[0]) {
			ctest.SomeVecs(test, indexes[0])
			ctempla = templa
		} else if templa.NVecs() > len(indexes[0]) {
			ctempla.SomeVecs(templa, indexes[0])
		} else {
			return nil, fmt.Errorf("chem.memMSD: Indexes don't match molecules")
		}
	} else {
		ctest.SomeVecs(test, indexes[0])
		ctempla.SomeVecs(templa, indexes[1])
	}
	if ctest.NVecs() != ctempla.NVecs() || tmp.NVecs() != ctest.NVecs() {
		return nil, fmt.Errorf("chem.memMSD: Ill formed matrices for memMSD calculation: cest: %d  ctempla: %d tmp: %d", ctest.NVecs(), ctempla.NVecs(), tmp.NVecs())
	}
	tmp.Sub(ctest, ctempla)
	msds := make([]float64, ctest.NVecs())
	for i := 0; i < tmp.NVecs(); i++ {
		v := tmp.VecView(i)
		msds[i] = v.Norm(2)
	}
	return msds, nil
}

//rMSD returns the RSMD (root of the mean square deviation) for the sets of cartesian
//coordinates in test and template, only considering the template and test atoms in
//the lists testlst and templalst, respectively. Since it is very explicit I leave it here for testing.
func rMSD(test, template *v3.Matrix, testlst, templalst []int) (float64, error) {
	//This is a VERY naive implementation.
	lists := [][]int{testlst, templalst}
	var ctest *v3.Matrix
	var ctempla *v3.Matrix
	if testlst == nil || len(testlst) == 0 {
		ctest = test
	} else {
		ctest = v3.Zeros(len(lists[0]))
		ctest.SomeVecs(test, lists[0])
	}
	if templalst == nil || len(templalst) == 0 {
		ctempla = template
	} else {
		ctempla = v3.Zeros(len(lists[1]))
		ctempla.SomeVecs(template, lists[1])
	}
	tmr, tmc := ctempla.Dims()
	tsr, tsc := ctest.Dims()
	if tmr != tsr || tmc != 3 || tsc != 3 {
		return -1, fmt.Errorf("RMSD: Ill formed matrices for RMSD calculation")
	}
	tr := tmr
	ctempla2 := v3.Zeros(ctempla.NVecs())
	ctempla2.Copy(ctempla)
	//the maybe thing might not be needed since we check the dimensions before.
	f := func() { ctempla2.Sub(ctempla2, ctest) }
	if err := gnMaybe(gnPanicker(f)); err != nil {
		return -1, CError{err.Error(), []string{"v3.Matrix.Sub", "RMSD"}}
	}
	var RMSD float64
	for i := 0; i < ctempla.NVecs(); i++ {
		temp := ctempla2.VecView(i)
		RMSD += math.Pow(temp.Norm(2), 2)
	}
	RMSD = RMSD / float64(tr)
	RMSD = math.Sqrt(RMSD)
	return RMSD, nil
}

//ImproperAlt calculates the improper dihedral between the points a, b,c,d
// as the angle between the plane defined by a,b,c and the cd vector
func ImproperAlt(a, b, c, d *v3.Matrix) float64 {
	all := []*v3.Matrix{a, b, c, d}
	for number, point := range all {
		pr, pc := point.Dims()
		if point == nil {
			panic(PanicMsg(fmt.Sprintf("goChem-Improper: Vector %d is nil", number)))
		}
		if pr != 1 || pc != 3 {
			panic(PanicMsg(fmt.Sprintf("goChem-Improper: Vector %d has invalid shape", number)))
		}
	}
	//bma=b minus a
	amb := v3.Zeros(1)
	cmb := v3.Zeros(1)
	dmc := v3.Zeros(1)
	amb.Sub(b, a)
	cmb.Sub(c, b)
	dmc.Sub(d, c)
	plane := cross(amb, cmb)
	angle := Angle(plane, dmc)
	return math.Pi/2.0 + angle
	/*	if angle <= math.Pi/2.0 {

		//	return math.Pi/2.0 - angle
			return math.Pi/2.0 + angle
		} else {
			//return (angle - math.Pi/2.0)
		}
	*/
}

//Improper calculates the improper dihedral between the points a, b,c,d
// as the angle between the plane defined by a,b,c and that defined by the plane bcd
func Improper(a, b, c, d *v3.Matrix) float64 {
	all := []*v3.Matrix{a, b, c, d}
	for number, point := range all {
		pr, pc := point.Dims()
		if point == nil {
			panic(PanicMsg(fmt.Sprintf("goChem-Improper: Vector %d is nil", number)))
		}
		if pr != 1 || pc != 3 {
			panic(PanicMsg(fmt.Sprintf("goChem-Improper: Vector %d has invalid shape", number)))
		}
	}
	//bma=b minus a
	amb := v3.Zeros(1)
	cmb := v3.Zeros(1)
	dmc := v3.Zeros(1)
	bmc := v3.Zeros(1)
	amb.Sub(a, b) //canged from Sub(b,a)
	cmb.Sub(c, b)
	bmc.Sub(b, c)
	dmc.Sub(d, c)
	plane1 := cross(amb, cmb)
	plane2 := cross(bmc, dmc)
	angle := Angle(plane1, plane2)
	return angle
	/*	if angle <= math.Pi/2.0 {

		//	return math.Pi/2.0 - angle
			return math.Pi/2.0 + angle
		} else {
			//return (angle - math.Pi/2.0)
		}
	*/
}

//DihedralAlt calculates the dihedral between the points a, b, c, d, where the first plane
//is defined by abc and the second by bcd.
//It is exactly the same as Dihedral, only kept for API stability.
func DihedralAlt(a, b, c, d *v3.Matrix) float64 {
	return Dihedral(a, b, c, d)
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
	amb := v3.Zeros(1)
	cmb := v3.Zeros(1)
	bmc := v3.Zeros(1)
	dmc := v3.Zeros(1)
	amb.Sub(a, b)
	cmb.Sub(c, b)
	bmc.Sub(b, c)
	dmc.Sub(d, c)
	v1 := cross(amb, cmb)
	v2 := cross(bmc, dmc)
	dihedral := Angle(v1, v2)
	//	dihedral := math.Atan2(first, second)
	return dihedral
}

//dihedral calculates the dihedral between the points a, b, c, d, where the first plane
//is defined by abc and the second by bcd.
//This is an old implemntation which does not give correct results
func dihedral(a, b, c, d *v3.Matrix) float64 {
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
	bmascaled.Scale(cmb.Norm(2), bma)
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
//linear and circular distortion, in that order, and error or nil.
//Based on the work of Taylor et al., .(1983), J Mol Graph, 1, 30
//This function has NOT been tested thoroughly in the sense of the appropiateness of the indexes definitions.
func RhoShapeIndexes(rhos []float64) (float64, float64, error) {
	if rhos == nil || len(rhos) < 3 {
		return -1, -1, CError{"goChe: Not enough or nil rhos", []string{"RhoShapeIndexes"}}
	}
	//	print(rhos[0],rhos[1],rhos[2]) ////////////////////////
	// Are these definitions reasonable?
	linear_distortion := (1 - (rhos[1] / rhos[0])) * 100     //Prolate
	circular_distortion := ((1 - (rhos[2] / rhos[1])) * 100) //Oblate
	return linear_distortion, circular_distortion, nil
}

//Rhos returns the semiaxis of the elipoid of inertia given the the moment of inertia tensor.
func Rhos(momentTensor *v3.Matrix, epsilon ...float64) ([]float64, error) {
	var e float64
	if len(epsilon) == 0 {
		e = -1
	} else {
		e = epsilon[0]
	}
	_, evals, err := v3.EigenWrap(momentTensor, e)
	if err != nil {
		return nil, errDecorate(err, "Rhos")
	}
	rhos := sort.Float64Slice{evals[0], evals[1], evals[2]} //invSqrt(evals[0]), invSqrt(evals[1]), invSqrt(evals[2])}
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
//if passed a nil Masser, it will simply set all masses to 1. If more than one Masser is passed
//Only the first will be considered
func BestPlane(coords *v3.Matrix, mol ...Masser) (*v3.Matrix, error) {
	var err error
	var Mmass []float64
	cr, _ := coords.Dims()
	if len(mol) != 0 && mol[0] != nil {
		Mmass, err = mol[0].Masses()
		if err != nil {
			return nil, errDecorate(err, "BestPlane")
		}
		if len(Mmass) != cr {
			return nil, CError{fmt.Sprintf("Inconsistent coordinates(%d)/atoms(%d)", len(Mmass), cr), []string{"BestPlane"}}
		}
	} else {
		Mmass = nil
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

//returns a float64 slice of the size requested filed with ones
func ones(size int) []float64 {
	slice := make([]float64, size, size)
	for k, _ := range slice {
		slice[k] = 1.0
	}
	return slice
}

//CenterOfMass returns the center of mass the atoms represented by the coordinates in geometry
//and the masses in mass, and an error. If no mass is given, it calculates the geometric center
func CenterOfMass(geometry *v3.Matrix, massS ...*mat.Dense) (*v3.Matrix, error) {
	var mass *mat.Dense
	if geometry == nil {
		return nil, CError{"goChem: nil matrix to get the center of mass", []string{"CenterOfMass"}}
	}
	gr, _ := geometry.Dims()
	if len(massS) == 0 || massS[0] == nil { //just obtain the geometric center
		tmp := ones(gr)
		mass = mat.NewDense(gr, 1, tmp) //gnOnes(gr, 1)
	} else {
		mass = massS[0]
	}
	tmp2 := ones(gr)
	gnOnesvector := mat.NewDense(1, gr, tmp2) //gnOnes(1, gr)

	ref := v3.Zeros(gr)
	ref.ScaleByCol(geometry, mass)
	ref2 := v3.Zeros(1)
	ref2.Mul(gnOnesvector, ref)
	ref2.Scale(1.0/mat.Sum(mass), ref2)
	return ref2, nil
}

//MassCenter centers in in the center of mass of oref. Mass must be
//A column vector. Returns the centered matrix and the displacement matrix.
func MassCenter(in, oref *v3.Matrix, massS ...*mat.Dense) (*v3.Matrix, *v3.Matrix, error) {
	or, _ := oref.Dims()
	ir, _ := in.Dims()
	var mass *mat.Dense
	if len(massS) == 0 || massS[0] == nil { //just obtain the geometric center
		tmp := ones(or)
		mass = mat.NewDense(or, 1, tmp) //gnOnes(or, 1)
	} else {
		mass = massS[0]
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
	ref2.Scale(1.0/mat.Sum(mass), ref2)
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
func MomentTensor(A *v3.Matrix, massslice ...[]float64) (*v3.Matrix, error) {
	ar, ac := A.Dims()
	var err error
	var mass *mat.Dense
	if len(massslice) == 0 || massslice[0] == nil {
		mass = gnOnes(ar, 1)
	} else {
		mass = mat.NewDense(ar, 1, massslice[0])
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
func AntiProjection(test, ref *v3.Matrix) *v3.Matrix {
	rr, _ := ref.Dims()
	testnorm := test.Norm(2)
	Uref := v3.Zeros(rr)
	Uref.Unit(ref)
	scalar := test.Dot(Uref)
	scalar = (testnorm * testnorm) / scalar
	Uref.Scale(scalar, Uref)
	return Uref
}
