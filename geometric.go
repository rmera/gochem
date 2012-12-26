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
import  "github.com/skelterjohn/go.matrix"
import "math"
import "sort"


//AngleInVectors takes 2 vectors and calculate the angle between them
//It does not check for correctness or return errors!
func AngleInVectors(v1,v2 *matrix.DenseMatrix) float64 {
	//Maybe I'll also write a safer version of this function?
	normproduct:=v1.TwoNorm()*v2.TwoNorm()
	dotprod:=Dot(v1,v2)
	argument:=dotprod/normproduct
	//Take care of floating point math errors
	if math.Abs(argument-1)<=appzero{
		argument=1
		}else if math.Abs(argument+1)<=appzero{
		argument=-1
		}
	//fmt.Println(dotprod/normproduct,argument) //dotprod/normproduct, dotprod, normproduct,v1.TwoNorm(),v2.TwoNorm())
	angle:=math.Acos(argument) 
	if math.Abs(angle)<=appzero{
		return 0.00
		}
	return angle
	}

/*
def angle_in_vectors(v1,v2): #calculates the angles between to vectors (Python Numeric arrays) in radians
	normproduct=norm(v1)*norm(v2)
	angle=np.arccos(np.dot(v1,v2)/normproduct)
	if angle<=approxzero or not angle:
		return 0
	return angle

 * */


//GetRotateToNewY takes a set of coordinates (mol) and a vector (y). It returns
//a rotation matrix that, when applied to mol, will rotate it around the Z axis 
//in such a way that the projection of newy in the XY plane will be aligned with
//the Y axis.
func GetRotateToNewY(mol, newy *matrix.DenseMatrix) (*matrix.DenseMatrix, error){
	if newy.Cols()!=3 || newy.Rows()!=1{
		return nil, fmt.Errorf("Wrong newy vector")
		}
	if mol.Cols()!=3{
		return nil, fmt.Errorf("Wrong mol vector")
		}
	gamma:=math.Atan2(newy.Get(0,0),newy.Get(0,1))
	singamma:=math.Sin(gamma)
	cosgamma:=math.Cos(gamma)
	operator:=[]float64{ cosgamma, singamma, 0,
	                    -singamma, cosgamma, 0,
	                        0,        0,    1}
	return matrix.MakeDenseMatrix(operator,3,3), nil
	                        
	
	}	


//GetRotateAroundZ takes a set of coordinates (mol) and a vector (y). It returns
//a rotation matrix that, when applied to mol, will rotate it around the Z axis 
//in such a way that the projection of newy in the XY plane will be aligned with
//the Y axis.
func GetRotateAroundZ(gamma float64) (*matrix.DenseMatrix, error){
	singamma:=math.Sin(gamma)
	cosgamma:=math.Cos(gamma)
	operator:=[]float64{ cosgamma, singamma, 0,
	                    -singamma, cosgamma, 0,
	                        0,        0,    1}
	return matrix.MakeDenseMatrix(operator,3,3), nil
	                        
	
	}	



//GetSwitchZ takes a matrix a row vector (newz).
//It returns a linear operator such that, when applied to a matrix mol ( with the operator on the right side)
//it will rotate mol such that the z axis is aligned with newz.
func GetSwitchZ(newz *matrix.DenseMatrix) (*matrix.DenseMatrix) {
	if newz.Cols()!=3 || newz.Rows()!=1{
		panic("Wrong newz vector")
		}
	normxy:=math.Sqrt(math.Pow(newz.Get(0,0),2)+math.Pow(newz.Get(0,1),2))
	theta:=math.Atan2(normxy,newz.Get(0,2)) //Around the new y
	phi:=math.Atan2(newz.Get(0,1),newz.Get(0,0))  //First around z
	psi:=0.000000000000  // second around z
	sinphi:=math.Sin(phi)
	cosphi:=math.Cos(phi)
	sintheta:=math.Sin(theta)
	costheta:=math.Cos(theta)
	sinpsi:=math.Sin(psi)
	cospsi:=math.Cos(psi)
	operator:=[]float64{cosphi*costheta*cospsi - sinphi*sinpsi,   -sinphi*cospsi - cosphi*costheta*sinpsi,   cosphi*sintheta,
				        sinphi*costheta*cospsi + cosphi*sinpsi,   -sinphi*costheta*sinpsi + cosphi*cospsi,   sintheta*sinphi,
				        -sintheta*cospsi,                          sintheta*sinpsi,                          costheta } 
	finalop:=matrix.MakeDenseMatrix(operator,3,3)
	return finalop
	
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
		return nil, nil, nil, nil, fmt.Errorf("GetSuper: Ill-formed matrices") 
		}
	var Scal float64
	p:=templa.Rows()
	Scal=float64(1.0)/float64(p)
	j:=matrix.Ones(p,1) //Mass is not important for this matter so we'll just use this.
	ctest,distest,err:=MassCentrate(test,test,j)
	if err!=nil{
		return nil, nil, nil, nil, err
		}
	ctempla,distempla,err:=MassCentrate(templa,templa,j)
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
	if err1 != nil {
		return nil, nil, nil, nil, err1
		}
	if  err2!=nil{
		return nil, nil, nil, nil, err2
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
	var RMSD float64
	for i:=0;i<ctempla.Rows();i++{
		temp:=ctempla.GetRowVector(i)
		RMSD+=math.Pow(temp.TwoNorm(),2)
		}
	RMSD=RMSD/float64(ctempla.Rows())
	RMSD=math.Sqrt(RMSD)
	return RMSD, nil
	}



//Dihedral calculate the dihedral between the points a, b, c, d, where the first plane 
//is defined by abc and the second by bcd.
func Dihedral(a,b,c,d *matrix.DenseMatrix) (float64){
	all:=[]*matrix.DenseMatrix{a,b,c,d}
	for number,point:=range(all){
		if point==nil{
			panic(fmt.Sprintf("Vector %d is nil",number))
			}
		if point.Rows()!=1 || point.Cols()!=3{
			panic(fmt.Sprintf("Vector %d has invalid shape",number))
			}
		}
	b1:=b.Copy()
	_=b1.Subtract(a)
	b2:=c.Copy()
	_=b2.Subtract(b)
	b3:=d.Copy()
	_=b3.Subtract(c)
	b1scaled:=b1.Copy()
	b1scaled.Scale(b2.TwoNorm())
	first:=Dot(b1scaled,Cross3DRow(b2,b3))
	second:=Dot(Cross3DRow(b1,b2),Cross3DRow(b2,b3))
	dihedral:=math.Atan2(first,second)                  
	return dihedral
	}

/***Shape indicator functions***/
const appzero float64 = 0.0000001  //used to correct floating point 
//errors. Everything equal or less than this is considered zero.


//point comparisons

//RhoShapeIndexes Get shape indices based on the axes of the elipsoid of inertia.
//Based on the work of Taylor et al., .(1983), J Mol Graph, 1, 30
//This function has NOT been tested.
func RhoShapeIndexes(evals []float64)(float64, float64, error){
	rhos,err:=Rhos(evals)
	linear_distortion:=(1-(rhos[1]/rhos[0]))*100 //Prolate
	circular_distortion:=(1-(rhos[2]/rhos[0]))*100 //Oblate
	return linear_distortion,circular_distortion,err
	}
	
//Rhos returns the semiaxis of the elipoid of inertia given the eigenvectors of the moment tensor.
func Rhos(evals []float64) ([]float64, error){
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


//BestPlaneP takes sorted evecs, according to the eval,s and returns a row vector that is normal to the
//Plane that best contains the molecule. Note that the function can't possibly check
//That the vectors are sorted!. The P at the end of the name is for Performance. If 
//That is not an issue it is safer to use the BestPlane function that wraps this one.
func BestPlaneP(evecs *matrix.DenseMatrix) (*matrix.DenseMatrix, error){
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


//BestPlane returns a row vector that is normal to the plane that best contains the molecule
//if passed a nil Ref, it will simply set all masses to 1.
func BestPlane(mol Ref, coords *matrix.DenseMatrix) (*matrix.DenseMatrix, error){
	var err error
	var Mmass *matrix.DenseMatrix
	if mol!=nil {
		if mol.Len()!=coords.Rows(){
			return nil, fmt.Errorf("Inconsistent coordinates(%d)/atoms(%d)",mol.Len(),coords.Rows())
			}
		Mmass,err=mol.MassCol()
		if err!=nil{
			return nil, err
			}
		}else{
		//Mmass=matrix.Ones(coords.Rows(),1)	
		}
	moment,err:=MomentTensor(coords,Mmass)
	if err!=nil{
		return nil, err
		}
	evecs,_,err:=Eigenwrap(moment)
	if err!=nil{
		return nil, err
		}
	normal,err:= BestPlaneP(evecs)
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
//	E.evecs[i],E.evecs[j]=E.evecs[j],E.evecs[i]
	E.evecs.SwapRows(i,j)
	}
func (E eigenpair)Len() int {
	return len(E.evals)
	}


//Eigenwrap wraps the matrix.DenseMatrix.Eigen() function in order to guarantee 
//That the eigenvectors and eigenvalues are sorted according to the eigenvalues
//It also guarantees orthonormality and handness. I don't know how many of 
//these are already guaranteed by Eig(). Will delete the unneeded parts 
//And even this whole function when sure.
func Eigenwrap(in *matrix.DenseMatrix) (*matrix.DenseMatrix, []float64, error){
	evecs,vals,_:=in.Eigen()
	evals:= [3]float64{vals.Get(0,0),vals.Get(1,1),vals.Get(2,2)}
	if err:=evecs.TransposeInPlace();err!=nil{
		return nil,nil,err
		}
	eig:=eigenpair{evecs,evals[:]}
	sort.Sort(eig)
	//Here I should orthonormalize vectors if needed instead of just complaining. 
	//I think orthonormality is guaranteed by  DenseMatrix.Eig() If it is, Ill delete all this
	//If not I'll add ortonormalization routines.
	for i:=0;i<eig.evecs.Rows();i++{
		vectori:=eig.evecs.GetRowVector(i)
		for j:=i+1;j<eig.evecs.Rows();j++{

	//		if i==j{
	//			continue //actually this cant happen, I could take this away.
	//			}
			if math.Abs(Dot(vectori,eig.evecs.GetRowVector(j)))>appzero{
				return eig.evecs,evals[:],fmt.Errorf("Vectors not ortogonal!")
				}
			}
		if math.Abs(vectori.TwoNorm()-1)>appzero{
			//Of course I could just normalize the vectors instead of complaining.
			return eig.evecs,evals[:],fmt.Errorf("Vectors not normalized")
			}
		}
	//Checking and fixing the handness of the matrix.This if-else is Jannes idea, 
	//I don't really know whether it works.
	eig.evecs.TransposeInPlace()
	if eig.evecs.Det()<0{
		eig.evecs.Scale(-1)
		} else {
		/*	
		eig.evecs.TransposeInPlace()
		eig.evecs.ScaleRow(0,-1)
		eig.evecs.ScaleRow(2,-1)
		eig.evecs.TransposeInPlace()
		*/
	//	fmt.Println("all good, I guess")
		}	
		eig.evecs.TransposeInPlace()
	return eig.evecs,eig.evals, nil  //Returns a slice of evals
	}



/*CenterOfMass returns the center of mass the atoms represented by the coordinates in geometry
and the masses in mass, and an error. If mass is nil, it calculates the geometric center*/
func CenterOfMass(geometry, mass *matrix.DenseMatrix)(*matrix.DenseMatrix,error){
	if geometry==nil{
		return nil, fmt.Errorf("nil matrix to get the center of mass")
		}

	if mass==nil{ //just obtain the geometric center
		mass=matrix.Ones(geometry.Rows(),1)
		}
	onesvector:=matrix.Ones(1,geometry.Rows())
	ref:=geometry.Copy()
	err:=DMScaleByCol(ref,mass)
	if err!=nil{
		return nil, err
		}
	ref2:=matrix.ParallelProduct(onesvector,ref)
	ref2.Scale(1.0/DMSummation(mass))
	return ref2,nil
	}


//MassCentrate centers in in the center of mass of oref. Mass must be
//A column vector. Returns the centered matrix and the displacement matrix.
func MassCentrate(in, oref, mass *matrix.DenseMatrix) (*matrix.DenseMatrix,*matrix.DenseMatrix, error){
	if mass==nil{ //just obtain the geometric center
		mass=matrix.Ones(oref.Rows(),1)
		}
	ref:=oref.Copy()
	onesvector:=matrix.Ones(1,ref.Rows())
	if err:=DMScaleByCol(ref,mass); err!=nil{
		return nil, nil, err
		}
	ref2:=matrix.ParallelProduct(onesvector,ref)
	ref2.Scale(1.0/DMSummation(mass))
	returned:=in.Copy()
	for i:=0;i<returned.Rows();i++{
		if err:=returned.GetRowVector(i).Subtract(ref2); err!=nil{
			return nil, nil, err
			}
		}
	return returned,ref2,nil
	}


//MomentTensor returns the moment tensor for a matrix A of coordinates and a column
//vector mass with the respective massess.
func MomentTensor(A, mass *matrix.DenseMatrix) (*matrix.DenseMatrix, error){
	if mass==nil{
		mass=matrix.Ones(A.Rows(),1)
		}
	center,_,err:=MassCentrate(A,A.Copy(),mass)
	if err!=nil{
		return nil, err
		}
	sqrmass:=DMPow(mass,0.5)
//	fmt.Println(center,sqrmass) ////////////////////////
	DMScaleByCol(center,sqrmass)
//	fmt.Println(center,"scaled center")
	centerT:=center.Transpose()
	moment:=matrix.ParallelProduct(centerT,center)
	return moment, err
	}


func Projection(test, ref *matrix.DenseMatrix) *matrix.DenseMatrix{
	Uref:=Unitarize(ref)
//	angle:=AngleInVectors(test,ref)
//	la:=test.TwoNorm()
	scalar:=Dot(test,Uref) //math.Abs(la)*math.Cos(angle)
	Uref.Scale(scalar)
	return Uref
	}

//Given a set of cartesian points in sellist, obtains a vector "plane" normal to the best plane passing through the points.
//It selects atoms from the set A that are inside a cone in the direction of "plane" that starts from the geometric center of the cartesian points,
//and has an angle of angle (radians), up to a distance distance. The cone is approximated by a set of radius-increasing cilinders with height thickness.
//If one starts from one given point, 2 cones, one in each direction, are possible. If whatcone is 0, both cones are considered.
//if whatcone<0, only the cone opposite to the plane vector direction. If whatcone>0, only the cone in the plane vector direction.
func SelCone(B, selection *matrix.DenseMatrix, angle, distance, thickness float64, whatcone int) []int{
	A:=B.Copy() //We will be altering the input so its better to work with a copy.
	selected:=make([]int,0,3)
	neverselected:=make([]int,0,300000) //waters that are too far to ever be selected
	nevercutoff:=distance/math.Cos(angle) //cutoff to be added to neverselected
	A,_,err:=MassCentrate(A,selection,nil)  //Centrate A in the geometric center of the selection, Its easier for the following calculations
	if err!=nil{
		panic(err.Error())
		}
	selection,_,_=MassCentrate(selection,selection,nil) //Centrate the selection as well
	plane,err:=BestPlane(nil,selection)  //I have NO idea which direction will this vector point. We might need its negative.
	if err!=nil{
		panic(err.Error())
	}
	for i:=thickness/2;i<=distance;i+=thickness{
		maxdist:=math.Tan(angle)*i  //this should give me the radius of the cone at this point
		for j:=0;j<A.Rows();j++{
			if isInInt(selected,j)  || isInInt(neverselected,j){ //we dont scan things that we have already selected, or are too far
				continue
			}
			atom:=A.GetRowVector(j)
			proj:=Projection(atom,plane)
			norm:=proj.TwoNorm()
			//Now at what side of the plane is the atom?
			angle:=AngleInVectors(atom,plane)
			if whatcone>0{
				if angle>math.Pi/2{
					continue
				}
			}else if whatcone<0{
				if angle<math.Pi/2{
					continue
				}	
			}
			if norm > i+(thickness/2.0) || norm < (i-thickness/2.0){
				continue
			}
			proj.Subtract(atom)
			projnorm:=proj.TwoNorm()
			if projnorm<=maxdist{
				selected=append(selected,j)
			}
			if projnorm>=nevercutoff{
				neverselected=append(neverselected,j)
			}
		}
	}
	return selected
}



