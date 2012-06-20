/*
 * geometric.go, part of gochem
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 * 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
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
import  "github.com/skelterjohn/go.matrix"
import "math"
import "sort"


	



//Geometrical functions.

//GetSwitchZ takes a matrix with cartesian coordinates in the rows (mol) and a row vector (newz).
//It returns a linear operator such that, when applied to the matrix mol (the operator on the right side)
//it will rotate mol such that the z axis is aligned with newz. Also returns error or nil.
//It appears to be wrong by 0.0014 radians in the X axis. Very little but I dont understand why. Could
//Just be numerical noise I guess.
func GetSwitchZ(mol, newz *matrix.DenseMatrix) (*matrix.DenseMatrix, error) {
	if newz.Cols()!=3 || newz.Rows()!=1{
		return nil, fmt.Errorf("Wrong newz vector")
		}
	if mol.Cols()!=3{
		return nil, fmt.Errorf("Wrong mol vector")
		}
	norm:=newz.TwoNorm()
	theta:=math.Atan2(norm,newz.Get(0,2)) //Around the new y
	phi:=math.Atan2(newz.Get(0,1),newz.Get(0,0))  //first around z
	psi:=0.000000000000  //Second turn around z
	sinphi:=math.Sin(phi)
	cosphi:=math.Cos(phi)
	sintheta:=math.Sin(theta)
	costheta:=math.Cos(theta)
	sinpsi:=math.Sin(psi)
	cospsi:=math.Cos(psi)
/*The operator used now is the transpose of this one
 * 	operator:=[]float64{cosphi*costheta*cospsi - sinphi*sinpsi,   sinphi*costheta*cospsi + cosphi*sinpsi, -sintheta*cospsi,
						-sinphi*cospsi - cosphi*costheta*sinpsi, -sinphi*costheta*sinpsi + cosphi*cospsi,  sintheta*sinpsi,
						cosphi*sintheta,                          sintheta*sinphi,                         costheta}
*/	
	//This one is the transpose of the other operator (commented). It is meant to be multiplied from the left by a row
	//vector.
	operator:=[]float64{cosphi*costheta*cospsi - sinphi*sinpsi,   -sinphi*cospsi - cosphi*costheta*sinpsi,   cosphi*sintheta,
				        sinphi*costheta*cospsi + cosphi*sinpsi,   -sinphi*costheta*sinpsi + cosphi*cospsi,   sintheta*sinphi,
				        -sintheta*cospsi,                          sintheta*sinpsi,                          costheta } 
	finalop:=matrix.MakeDenseMatrix(operator,3,3)
	fmt.Println(operator, "\n\n",finalop)
	return finalop, nil
	
	}




//Super determines the best rotation and translations to superimpose the atoms of molecule test, frame frametest,
//listed in testlst on te atoms of molecule templa, frame frametempla, listed in templalst. 
//It applies those rotation and translations to the whole frame frametest of molecule test, in palce. 
//testlst and templalst must have the same number of elements.
func Super(test, templa *Molecule, testlst, templalst []int, frametest, frametempla int) error{
	ctest,err1:=test.GetCoords(testlst,frametest)
	ctempla,err2:=templa.GetCoords(templalst,frametempla)
	if err1!=nil || err2!=nil{
		return fmt.Errorf("Frame numbers given for test or template out of range")
		}
	_,rotation,trans1,trans2,err1:=GetSuper(ctest,ctempla)
	if err1!=nil{
		return err1
		}
	err1=AddRow(test.Coords[frametest],trans1)
	test.Coords[frametest]=matrix.ParallelProduct(test.Coords[frametest],rotation)
	err2=AddRow(test.Coords[frametest],trans2)
	if err1 != nil || err2!=nil{
		return fmt.Errorf("Unexpected error when aplying superposition")
		}
	return nil
	}


//GetSuper superimposes the set of cartesian coordinates given as the rows of the matrix test on the ones of the rows
//of the matrix templa. Returns the transformed matrix, the rotation matrix, 2 translation row vectors
//For the superposition plus an error. In order to perform the superposition, without using the transformed
//the first translation vector has to be added first to the moving matrix, then the rotation must be performed
//and finally the second translation has to be added.
//This is a low level function, although one can use it directly since it returns the transformed matrix.
//The math for this function is by Prof. Veronica Jimenez-Curihual, University of Concepcion, Chile.
func GetSuper(test, templa *matrix.DenseMatrix)(*matrix.DenseMatrix, *matrix.DenseMatrix, *matrix.DenseMatrix, *matrix.DenseMatrix, error){
	dot:=matrix.ParallelProduct
	if templa.Rows() != test.Rows() || templa.Cols()!= 3 || test.Cols()!=3{
		return nil, nil, nil, nil, fmt.Errorf("Ill-formed matrices") 
		}
	var Scal float64
	p:=templa.Rows()
	Scal=float64(1.0)/float64(p)
	j:=matrix.Ones(p,1) //Mass is not important for this matter so we'll just use this.
	ctest,distest,err:=CenterMass(test,test,j)
	if err!=nil{
		return nil, nil, nil, nil, err
		}
	ctempla,distempla,err:=CenterMass(templa,templa,j)
	if err!=nil{
		return nil, nil, nil, nil, err
		}
	Mid:=matrix.Eye(p)
	jT:=j.Transpose()
	ScaledjProd:=dot(j,jT) 
	ScaledjProd.Scale(Scal)
	Maux:=dot(dot(ctempla.Transpose(),Mid),ctest)
	Maux=Maux.Transpose() //Dont understand why this is needed
	U,_,Vt,err:=Maux.SVD()
	if err!=nil{
		return nil, nil, nil, nil, err
		}
	U.Scale(-1)
	Vt.Scale(-1)
	//SVD gives different results here than in numpy. U and Vt are multiplide by -1 in one of them
	//and gomatrix gives as Vt the transpose of the matrix given as Vt by numpy. I guess is not an 
	//error, but don't know for sure.
	Rotation:=dot(Vt,U.Transpose())
	Rotation=Rotation.Transpose() //Don't know why does this work :(
	if Rotation.Det()<0{
		return nil, nil, nil, nil,fmt.Errorf("Got a reflection instead of a translations. The objects may be specular images of each others")
		}
	jT.Scale(Scal)
	subtempla:=matrix.MakeDenseCopy(ctempla)
	subtempla.SubtractDense(dot(ctest,Rotation))
	Translation:=dot(jT,subtempla)
	err1:=Translation.Add(distempla)
	//This allings the transformed with the original template, not the mean centrate one
	transformed:=matrix.Product(ctest,Rotation)
	err2:=AddRow(transformed,Translation)
	//end transformed
	distest.Scale(-1)
	if err1 != nil || err2!=nil{
		return nil, nil, nil, nil, fmt.Errorf("Problem with the final translations in superposition procedure")
		}
	return transformed, Rotation, distest, Translation, nil
	}

//I keep this just in case I manage to fix it at some point
func rmsd_fail(test, template *matrix.DenseMatrix) (float64, error){
	ctempla:=template.Copy()
	err:=ctempla.Subtract(test)
	if err!=nil{
		return 9999999, err
		}
	dev:=matrix.ParallelProduct(ctempla.Transpose(),ctempla)
	RMSDv:=dev.Trace()
	RMSDv=math.Sqrt(RMSDv)
	return RMSDv, nil
	}

//RMSD returns the RSMD (root of the mean square deviation) for the sets of cartesian
//coordinates in test and template
func RMSD(test, template *matrix.DenseMatrix) (float64, error){
	if template.Rows()!=test.Rows() || template.Cols()!= 3 || test.Cols()!=3{
		return 0, fmt.Errorf("Ill formed matrices for RMSD calculation")
		}
	ctempla:=template.Copy()
	err:=ctempla.Subtract(test)
	if err!=nil{
		return 0, err
		}	
	DMPowInPlace(ctempla,2)
	var RMSD float64
	for i:=0;i<ctempla.Rows();i++{
		temp:=ctempla.GetRowVector(i)
		RMSD+=temp.TwoNorm()
		}
	RMSD=math.Sqrt(RMSD)
	return RMSD, nil
	}



/***Shape indicator functions***/
const appzero float64 = 0.000001  //used at least in Eigenwrap to make floating
//point comparisons

//GetRhoShapeIndexes Get shape indices based on the axes of the elipsoid of inertia.
//Based on the work of Taylor et al., .(1983), J Mol Graph, 1, 30
//This function has NOT been tested.
func GetRhoShapeIndexes(evals []float64)(float64, float64, error){
	rhos,err:=GetRhos(evals)
	linear_distortion:=(1-(rhos[1]/rhos[0]))*100 //Prolate
	circular_distortion:=(1-(rhos[2]/rhos[0]))*100 //Oblate
	return linear_distortion,circular_distortion,err
	}
	
//GetRhos returns the semiaxis of the elipoid of inertia given the eigenvectors of the moment tensor.
func GetRhos(evals []float64) ([]float64, error){
	rhos:=sort.Float64Slice{InvSqrt(evals[0]),InvSqrt(evals[1]),InvSqrt(evals[2])}
	if evals[2]<=appzero{
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


//GetBestPlaneB takes sorted evecs, according to the eval,s and returns a row vector that is normal to the
//Plane that best contains the molecule. Note that the function can't possibly check
//That the vectors are sorted!. The P at the end of the name is for Performance. If 
//That is not an issue it is safer to use the GetBestPlane function that wraps this one.
func GetBestPlaneP(evecs *matrix.DenseMatrix) (*matrix.DenseMatrix, error){
	if evecs.Rows()!=3 || evecs.Cols()!=3{
		return evecs, fmt.Errorf("Eigenvectors matrix must be 3x3")
		}
	v1:=evecs.GetColVector(2)
	v2:=evecs.GetColVector(1)
	tv1:=v1.Transpose()
	tv2:=v2.Transpose()
	normal:=Cross3DRow(tv1,tv2)
	return normal, nil
	}

//GetBestPlane returns a row vector that is normal to the plane that best contains the molecule
func GetBestPlane(mol *Molecule, frame int,masses bool) (*matrix.DenseMatrix, error){
	var Mmass *matrix.DenseMatrix
	if len(mol.Atoms)!=mol.Coords[frame].Rows(){
		return nil, fmt.Errorf("Inconsistent coordinates/atoms in frame %d", frame)
		}
	if masses {
		mass,err:=mol.GetMassArray()
		if err!=nil{
			return nil, err
			}
		Mmass=matrix.MakeDenseMatrix(mass[:],1,len(mass))
		}else{
		Mmass=matrix.Ones(1,len(mol.Atoms))	
		}
	moment,err:=MomentTensor(mol.Coords[frame],Mmass)
	if err!=nil{
		return nil, err
		}
	evecs,_,err:=Eigenwrap(moment)
	if err!=nil{
		return nil, err
		}
	normal,err:= GetBestPlaneP(evecs)
	//MomentTensor(, mass) 
	return normal, err
	}
	 




//This is a facility to sort Eigenvectors/Eigenvalues pairs
//It satisfies the sort.Interface interface.
type eigenpair struct{
	//evecs must have as many rows as evals has elements.
	evecs *matrix.DenseMatrix 
	evals sort.Float64Slice
	}
func (E eigenpair)Less(i,j int) bool {
	return E.evals[i]<E.evals[j]
	}
func (E eigenpair)Swap(i,j int) {
	E.evals.Swap(i,j)
	E.evecs.SwapRows(i,j)
	}
func (E eigenpair)Len() int {
	return len(E.evals)
	}


//Eigenapir wraps the matrix.DenseMatrix.Eigen() function in order to guarantee 
//That the eigenvectors and eigenvalues are sorted according to the eigenvalues
//and also orthonormality and Handness I don't know how many of these are already 
//guaranteed by Eig(). Will delete the unneeded parts when sure.
func Eigenwrap(in *matrix.DenseMatrix) (*matrix.DenseMatrix, []float64, error){
	evecs,vals,_:=in.Eigen()
//	evecs:= [3]*matrix.DenseMatrix{vecs.GetRowVector(0),vecs.GetRowVector(1),vecs.GetRowVector(2)}
	evals:= [3]float64{vals.Get(0,0),vals.Get(1,1),vals.Get(2,2)}
	eig:=eigenpair{evecs,evals[:]}
	sort.Sort(eig)
	//Here I should orthonormalize vectors if needed instead of just complaining. 
	//I think orthonormality is guaranteed by  DenseMatrix.Eig() If it is, Ill delete all this
	//If not I'll add ortonormalization routines.
	for i:=0;i<eig.evecs.Rows();i++{
		vectori:=eig.evecs.GetRowVector(i)
		for j:=i+1;j<eig.evecs.Rows();j++{
			if i==j{
				continue
				}
			if math.Abs(Dot(vectori,eig.evecs.GetRowVector(j)))>appzero{
				return eig.evecs,evals[:],fmt.Errorf("Vectors not ortogonal!")
				}
			}
		if math.Abs(vectori.TwoNorm()-1)>appzero{
			return eig.evecs,evals[:],fmt.Errorf("Vectors not normalized")
			}
		}
	//Checking and fixing the handness of the matrix.This if-else is Jannes idea, 
	//I don't really know whether it works.
	if eig.evecs.Det()<0{
		eig.evecs.Scale(-1)
		} else {
		eig.evecs.TransposeInPlace()
		eig.evecs.ScaleRow(0,-1)
		eig.evecs.ScaleRow(2,-1)
		eig.evecs.TransposeInPlace()
		}	
	return eig.evecs,evals[:], nil  //Returns a slice of evals
	}




//CenterMass centers in in the center of mass of ref. Mass must be
//A column vector. Returns the centered matrix and the displacement matrix.
func CenterMass(in, oref, mass *matrix.DenseMatrix) (*matrix.DenseMatrix,*matrix.DenseMatrix, error){
	ref:=oref.Copy()
	onesvector:=matrix.Ones(1,in.Rows())
	err:=DMScaleByCol(ref,mass)
	if err!=nil{
		return nil, nil, err
		}
	ref2:=matrix.ParallelProduct(onesvector,ref)
	ref2.Scale(1.0/DMSummation(mass))
	returned:=in.Copy()
	for i:=0;i<returned.Rows();i++{
		err=returned.GetRowVector(i).Subtract(ref2)
		if err!=nil{
			return nil, nil, err
			}
		}
	return returned,ref2,err
	}


//MomentTensor returns the moment tensor for a matrix A of coordinates and a column
//vector mass with the respective massess.
func MomentTensor(A, mass *matrix.DenseMatrix) (*matrix.DenseMatrix, error){
	center,_,err:=CenterMass(A,A.Copy(),mass)
	if err!=nil{
		return nil, err
		}
	sqrmass:=DMPow(mass,0.5)
	DMScaleByCol(center,sqrmass)
	moment:=matrix.ParallelProduct(center.Transpose(),center)
	return moment, err
	}



//Theses ones are basic math, belongs more to the go.matrix package
//If there is something similar already made
//in go.matrix this functions will be deleted. Otherwise they could be
//made methods for DenseMatrix and included in go.matrix

//Some of this functions don't return error messages because they are meant to
//Be inserted in mathematical expressions and thus they need to return only one value.


//AddRow adds the row vector row to each row of the matrix big, in place. Both need the same ammount of columns.
func AddRow(big,row *matrix.DenseMatrix)(error){
	bigrows:=big.Rows()
	if big.Cols() != row.Cols() || row.Rows()!=1{
		return fmt.Errorf("Ill-formed matrices for multiplication")
		}
	for i:=0;i<bigrows;i++{
		j:=big.GetRowVector(i)
		j.Add(row)
		}
	return nil
	}

//Cross3D Takes 2 3-len column or row vectors and returns a column or a row
//vector, respectively, with the Cross product of them.
func Cross3D(a,b *matrix.DenseMatrix)(*matrix.DenseMatrix,error){
	ac:=a.Cols()
	ar:=a.Rows()
	bc:=b.Cols()
	br:=b.Rows()
	if ac != bc || ar != br {
		return nil, fmt.Errorf("ill-formed vectors for cross product")
		}
	if ac!=3 {
		//Ok, Im sure one can do this better.
		c:=a.Transpose()
		d:=b.Transpose()
		e:=Cross3DRow(c,d)
		f:=e.Transpose()
		return f, nil
		}
	if ar!=3 {
		return nil, fmt.Errorf("Malformed vectors for cross product")
		}
	return Cross3DRow(a,b), nil
	}

//Cross3DRow returns the cross product of 2 row vectors. No error checking!
func Cross3DRow(a,b *matrix.DenseMatrix)*matrix.DenseMatrix{
	vec:=make([]float64,3,3)
	vec[0]=a.Get(0,1)*b.Get(0,2) - a.Get(0,2)*b.Get(0,1)
	vec[1]=a.Get(0,2)*b.Get(0,0) - a.Get(0,0)*b.Get(0,2)
	vec[1]=a.Get(0,0)*b.Get(0,1) - a.Get(0,1)*b.Get(0,0)
	return matrix.MakeDenseMatrix(vec,1,3)
	}


//InvSqrt return the inverse of the square root of val, or zero if
//val<appzero. It doesn't check for negative numbers! 
func InvSqrt(val float64) float64{
	if val<=appzero{
		return 0
		}
	return 1.0/math.Sqrt(val)	
	}

//KronekerDelta is a naive implementation of the kroneker delta function.
func KronekerDelta(a,b float64) float64{
	if math.Abs(a-b)<=appzero{
		return 1
		}
	return 0	
	}
	
	

//Dot returns the dot product between 2 vectors or matrices. Just the sum of the 
//Element-wise multiplication. In this case returning error
//makes it problematic to use in complex operations, so this returns -1
//when problems.
func Dot(A, B *matrix.DenseMatrix) float64{
	var err error
	if A.Cols()!=B.Cols() || A.Rows()!=B.Rows(){
		return -1
		}
	//For some crazy reason if the F variable is called C, I get a
	//"ScaleMatrixDense undeclared" error at compile time :S
	F := A.Copy()
	err=F.ScaleMatrixDense(B)
	if err!=nil{
		return -1
		}
	return DMSummation(F)
	}

//DMPowInPlace raises the DenseMatrix A, element-wise, to the nth power.  
func DMPowInPlace(A *matrix.DenseMatrix, n float64){
	for i:=0;i<A.Rows();i++{
		for j:=0;j<A.Cols();j++{
		//	fmt.Println(i,j, A.Rows(),A.Cols())
			A.Set(i,j,math.Pow(A.Get(i,j),n))
			}
		}
	}

//DMPow returns A^n. It does not modify A.
func DMPow(B *matrix.DenseMatrix, n float64) *matrix.DenseMatrix{
	A:=B.Copy()
	for i:=0;i<A.Rows();i++{
		for j:=0;j<A.Cols();j++{
		//	fmt.Println(i,j, A.Rows(),A.Cols())
			A.Set(i,j,math.Pow(A.Get(i,j),n))
			}
		}
	return B
	}

//DMScaleByCol scales each column of matrix A by Col.
func DMScaleByCol(A, Col *matrix.DenseMatrix) error{
	Rows:=A.Rows()
	if Rows != Col.Rows() || Col.Cols()>1{
		return fmt.Errorf("Malformed matrices for scaling")
		}
	for i:=0;i<Rows;i++{
		A.ScaleRow(i,Col.Get(i,0))
		}
	return nil
	}

//DMScaleByRow each row of matrix A by Row. Only for square matrices.
func DMScaleByRow(A, Row *matrix.DenseMatrix) error{
	A.TransposeInPlace()  //I guess this is not SUPER efficient
	Col:=Row.Transpose()
	err:=DMScaleByCol(A, Col)
	A.TransposeInPlace()
	return err
	}

//DMSummation returns the sum of all elements in matrix A.
func DMSummation(A *matrix.DenseMatrix) float64 {
	Rows:=A.Rows()
	Cols:=A.Cols()
	var sum float64
	for i:=0;i<Cols;i++{
		for j:=0;j<Rows;j++{
			sum+=A.Get(j,i)
			}
		}
	return sum
	}
	
//DMaddfloat returns a matrix which elements are those of matrix A plus the float B.
func DMaddfloat(A *matrix.DenseMatrix, B float64)*matrix.DenseMatrix{
	Rows:=A.Rows()
	Cols:=A.Cols()
	copy:=matrix.MakeDenseCopy(A)
	for i:=0;i<Cols;i++{
		for j:=0;j<Rows;j++{
			copy.Set(j,i,(A.Get(j,i)+B))
			}
		}
	return copy
	}






