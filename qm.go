/*
 * qm.go, part of gochem.
 * 
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
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
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.  
 * 
 * 
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

	
package chem

import "os"  
import "strings"
import "fmt"
import  "github.com/skelterjohn/go.matrix"

type IntConstraint struct{
	Kind byte
	Atoms []int
	}

type  PointCharge struct{
	Charge float64
	Coords *matrix.DenseMatrix
	}

type QMCalc struct {
	Method string
	Basis string
	AuxBasis string //for RI calculations
	AuxColBasis string //for RICOSX or similar calculations
	HighBasis string  //a bigger basis for certain atoms
	Lowbasis string  //lower basis for certain atoms
	HBatoms []int
	LBatoms []int
	HBelements []string
	LBelements []string
	CConstraint []int //cartesian contraints
	IConstraints []IntCostraint //internal constraints
	Dielectric float64
	Solventmethod string
	Disperssion string  //D, D2, D3
	Others string //analysis methods, etc
	PCharges []PointCharge
	Guess string //initial guess
	OldMO bool //Try to look for a file with MO. The 
	
	}


type OrcaRunner struct{
	defmethod string
	defbasis string
	tightness string
	slowconv string
	previousMO string
	command string
	inputname string
	nCPU int
	optimize bool
	}

func MakeOrcaRunner() *OrcaRunner{
	run:=new(OrcaRunner)
	run.SetDefaults()
	return run 
	}

func (O *OrcaRunner)SetDefaults(){
	O.defmethod="revPBE"
	O.defbasis="def2-SVP"
	O.defauxbasis="def2-SVP/J"
	O.command=os.ExpandEnv("${ORCA_PATH}/orca")
	if O.command=="/orca"{ //if ORCA_PATH was not defined
		O.command="./orca"
		}
	}

//BuildInput builds an input for ORCA based int the data in atoms, coords and C.
//returns only error.
func (O *OrcaRunner) BuildInput(atoms Reference, coords *matrix.DenseMatrix, C *QMCalc) error{
	//set defaults
	if Q.Basis==""{
		fmt.Fprintf(os.Stderr,"no basis set assigned for ORCA calculation, will used the default %s, \n",O.defmethod)
		Q.Basis=O.defbasis	
	if Q.Method==""{
		fmt.Fprintf(os.Stderr,"no method assigned for ORCA calculation, will used the default %s, \n",O.defmethod)
		Q.Method=O.defmethod
		Q.AuxColBasis="" //makes no sense for pure functional
		Q.AuxBasis=Sprintf("%s/J",Q.Basis)
		}
	//The usage of RI/RICOSX is given by the presence of AuxBasis/AuxColBasis. Its the user responsability 
	//not to set AuxColBasis for non-hybrid functionals and not to set AuxBasis and not AuxColBasis for
	//hybrid functionals.
	if Q.AuxColBasis!=""{
		//this will fail if someone wrote RI in the Other variable.
		if !(strings.Contains("RICOSX",Q.Others)){
			Q.Others=Sprintf("%s %s",Q.Others,"RICOSX")
			}
		}else if Q.AuxBasis!=""{
		if !(strings.Contains("RI",Q.Others)){
			Q.Others=Sprintf("%s %s",Q.Others,"RI")
			}
		}
	if Q.Method=="" || basis==""{
		return fmt.Errorf("Not enough options for the optimization")
		}
	disp=""
	switch Q.Disperssion{
		case "":
		case "D2":
		disp="VDW"
		case "D3":
		disp="VDW3"
		default:
		disp="VDW3"
		}
	if atoms==nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
		}
	opt:=""
	if O.optimise==true{
		opt=="Opt"
		}
	//If this flag is set we'll look for a suitable MO file.
	//If not found, we'll just use the default ORCA guess
	hfuhf="RHF"
	if atoms.Unpaired()!=0{
		hfufh="UHF"
		}
	if O.OldMO==true{
		files:=os.ReadDir(-1) //Get all the files
		moname:=O.previousMO
		O.previousMO=""
		for _,val:= range(files){
			if val.IsDir()==true{
				continue
				}
			name:=val.Name()
			if moname!=""{
				if name==moname{
					O.previousMO=name
					break
					}
				}else if string.Contains(".gbw",name){
				O.previousMO=name
				break
				}
			}
		if O.previousMO!=""{	
			Q.Guess="MORead"
			moinp:=fmt.Sprintf('\%moinp "%s\n"',O.previousMO)
			}else{
			moinp:=""
			Q.Guess=""	//The default guess
			}
		}
	if O.tightness==""{
		tight="TightSCF"
		}else{
		tight=O.tightness	
		}
	if slowconv==""{
		//The default for this is nothing for RHF and SlowConv for UHF
		if atoms.Unpaired>0{
			slowconv="SlowConv"
			}
		}
	pal:=""
	if O.nCPU>1{
		if O.nCPU>8{
			return fmt.Errorf("CPU number of %d for ORCA calculations currently not supported, maximun 8", O.nCPU)
			}
		pal:=fmt.Sprintf("PAL%d",O.nCPU)
		}
	MainOptions:=[]string{"!",hfuhf,Q.Method,Q.basis,Q.AuxBasis,Q.AuxColBasis,tight,disp,Q.Guess,opt,Q.Others,pal,"\n"}
	mainline:=strings.Join(MainOptions," ")
	constraints:=O.buildCCconstraints(Q.CConstraints)
	cosmo=""
	if C.Dielectric>=0{
		cosmo:=fmt.Sprintf("\%cosmo epsilon %d\n        refrac 1.30\n        end\n",C.Dielectric)
		}
	ElementBasis:=""
	if C.HBelements!=nil || C.LBelements!=nil{
		elementbasis:=make([]string,0,len(C.HBelements)+len(C.LBelements)+2)
		elementbasis:=append(elementbasis," \%basis \n")		
		for _,val:=range(C.HBelements){
			elementbasis:=append(elementbasis,fmt.Sprintf('  newgto %s "%s" end\n',val,HighBasis))
			}
		for _,val:= range(C.LBelements){
			elementbasis:=append(elementbasis,fmt.Sprintf('  newgto %s "%s" end\n',val,LowBasis))
			}
		elementbasis:=append(elementbasis,"         end\n")
		ElementBasis=strings.Join(elementbasis,"")
		}
	//Now lets write the thing
	if O.inputname==nil{
		O.inputname="Gochem.inp"
		}
	file,err:=os.Create(O.inputname)
	if err!=nil{
		return err
		}
	defer file.Close()
	_,err:=fmt.Fprint(file,mainline)
	//With this check its assumed that the file is ok.
	if err!=nil{
		return err
		}
	fmt.Fprint(file,moinp)
	fmt.Fprint(file,constraints)
	fmt.Fprint(file,ElementBasis)
	fmt.Fprint(file,cosmo)
	fmt.Fprint(file,"\n")
	//Now the type of coords, charge and multiplicity
	fmt.Fprint(file,"* xyz %d %d\n",atoms.Charge(),atoms.Unpaired()+1)
	//now the coordinates
	for i:=i<atoms.Len();i++{
		fmt.Fprintf(file,"%-2s  %8.3f%8.3f%8.3f \n",mol.Atom(i).Symbol, coords.Get(i,0), coords.Get(i,2), coords.Get(i,2))	
		}
	fmt.Fprintf(file,"*\n")
	
	}
	
//buildCConstraints transforms the list of cartesian constrains in the QMCalc structre
//into a string with ORCA-formatted cartesian constraints
(O *OrcaRunner) buildCConstraints(C []int) string{
	if C==nil{
		return "\n" //no constraints
		}
	constraints:=make([]string,len(C)+3)
	constraints[0]="\%geom Constraints\n"
	for key,val:=C{
		constraints[key+1]=fmt.Sprintf("         {C %d C}\n",val)
		}
	last:=len(constraints)-1
	constraints[last-1]="         end\n"
	constraints[last]=" end\n"
	final:=strings.join(constraints,"")
	return final
	}







