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
import "runtime"
import "os/exec"

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
	LowBasis string  //lower basis for certain atoms
	HBAtoms []int
	LBAtoms []int
	HBElements []string
	LBElements []string
	CConstraints []int //cartesian contraints
//	IConstraints []IntConstraint //internal constraints
	Dielectric float64
//	Solventmethod string
	Disperssion string  //D, D2, D3
	Others string //analysis methods, etc
//	PCharges []PointCharge
	Guess string //initial guess
	OldMO bool //Try to look for a file with MO. The 
	Optimize bool
	SCFTightness int
	SCFConvHelp int
	}


type OrcaRunner struct{
	defmethod string
	defbasis string
	defauxbasis string
	previousMO string
	command string
	inputname string
	nCPU int
	}



func MakeOrcaRunner() *OrcaRunner{
	run:=new(OrcaRunner)
	run.SetDefaults()
	return run 
	}


func (O *OrcaRunner)SetnCPU(cpu int){
	O.nCPU=cpu
	}


func (O *OrcaRunner)SetInputName(name string){
	O.inputname=name
	}

func (O *OrcaRunner)SetCommand(name string){
	O.command=name
	}



func (O *OrcaRunner)SetDefaults(){
	O.defmethod="revPBE"
	O.defbasis="def2-SVP"
	O.defauxbasis="def2-SVP/J"
	O.command=os.ExpandEnv("${ORCA_PATH}/orca")
	if O.command=="/orca"{ //if ORCA_PATH was not defined
		O.command="./orca"
		}
	cpu:=runtime.NumCPU()
	if cpu>8{
		O.nCPU=8
		}
	O.nCPU=cpu
	
	}

//BuildInput builds an input for ORCA based int the data in atoms, coords and C.
//returns only error.
func (O *OrcaRunner) BuildInput(atoms Ref, coords *matrix.DenseMatrix, Q *QMCalc) error{
	//Only error so far
	if atoms==nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
		}
	if Q.Basis==""{
		fmt.Fprintf(os.Stderr,"no basis set assigned for ORCA calculation, will used the default %s, \n",O.defbasis)
		Q.Basis=O.defbasis	
		}
	if Q.Method==""{
		fmt.Fprintf(os.Stderr,"no method assigned for ORCA calculation, will used the default %s, \n",O.defmethod)
		Q.Method=O.defmethod
		Q.AuxColBasis="" //makes no sense for pure functional
		Q.AuxBasis=fmt.Sprintf("%s/J",Q.Basis)
		}
	//The usage of RI/RICOSX is given by the presence of AuxBasis/AuxColBasis. Its the user responsability 
	//not to set AuxColBasis for non-hybrid functionals and not to set AuxBasis and not AuxColBasis for
	//hybrid functionals.
	if Q.AuxColBasis!=""{
		//this will fail if someone wrote RI in the Other variable.
		if !(strings.Contains("RICOSX",Q.Others)){
			Q.Others=fmt.Sprintf("%s %s",Q.Others,"RICOSX")
			}
	
		}else if Q.AuxBasis!=""{
			Q.Others=fmt.Sprintf("%s %s",Q.Others,"RI")
			}
	disp:="VDW3"
	if Q.Disperssion!=""{
		disp=orcaDisp[Q.Disperssion]	
		}
	opt:=""
	if Q.Optimize==true{
		opt="Opt"
		}
	//If this flag is set we'll look for a suitable MO file.
	//If not found, we'll just use the default ORCA guess
	hfuhf:="RHF"
	if atoms.Unpaired()!=0{
		hfuhf="UHF"
		}
	moinp:=""
	if Q.OldMO==true{
		dir,_:=os.Open("./") //This should always work, hence ignoring the error
		files,_:=dir.Readdir(-1) //Get all the files. 
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
				}else if strings.Contains(".gbw",name){
				O.previousMO=name
				break
				}
			}
		if O.previousMO!=""{	
			Q.Guess="MORead"
			moinp=fmt.Sprintf("%moinp \"%s\n\"",O.previousMO)
			}else{
			moinp=""
			Q.Guess=""	//The default guess
			}
		}
	tight:="TightSCF"
	if Q.SCFTightness!=0{
		tight=orcaSCFTight[Q.SCFTightness]
		}
	conv:=""
	if Q.SCFConvHelp==0{
		//The default for this is nothing for RHF and SlowConv for UHF
		if atoms.Unpaired()>0{
			conv="SlowConv"
			}
		}else{
		conv=orcaSCFConv[Q.SCFConvHelp]
		}
	pal:=""
	if O.nCPU>1{
		if O.nCPU>8{
			fmt.Fprintf(os.Stderr,"CPU number of %d for ORCA calculations currently not supported, maximun 8", O.nCPU)
			O.nCPU=8
			}
		pal=fmt.Sprintf("PAL%d",O.nCPU)
		}
	MainOptions:=[]string{"!",hfuhf,Q.Method,Q.Basis,Q.AuxBasis,Q.AuxColBasis,tight,disp,conv,Q.Guess,opt,Q.Others,pal,"\n"}
	mainline:=strings.Join(MainOptions," ")
	constraints:=O.buildCConstraints(Q.CConstraints)
	cosmo:=""
	if Q.Dielectric>0{
		cosmo=fmt.Sprintf("%%cosmo epsilon %1.0f\n        refrac 1.30\n        end\n",Q.Dielectric)
		}
	ElementBasis:=""
	if Q.HBElements!=nil || Q.LBElements!=nil{
		elementbasis:=make([]string,0,len(Q.HBElements)+len(Q.LBElements)+2)
		elementbasis=append(elementbasis," %basis \n")		
		for _,val:=range(Q.HBElements){
			elementbasis=append(elementbasis,fmt.Sprintf("  newgto %s \"%s\" end\n",val,Q.HighBasis))
			}
		for _,val:= range(Q.LBElements){
			elementbasis=append(elementbasis,fmt.Sprintf("  newgto %s \"%s\" end\n",val,Q.LowBasis))
			}
		elementbasis=append(elementbasis,"         end\n")
		ElementBasis=strings.Join(elementbasis,"")
		}
	//Now lets write the thing
	if O.inputname==""{
		O.inputname="Gochem"
		}
	
	file,err:=os.Create(fmt.Sprintf("%s.inp",O.inputname))
	if err!=nil{
		return err
		}
	defer file.Close()
	_,err=fmt.Fprint(file,mainline)
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
	fmt.Fprintf(file,"* xyz %d %d\n",atoms.Charge(),atoms.Unpaired()+1)
	//now the coordinates
	for i:=0;i<atoms.Len();i++{
		newbasis:=""
		if isInInt(Q.HBAtoms,i)==true{
			newbasis=fmt.Sprintf("newgto \"%s\" end",Q.HighBasis)
			}else if isInInt(Q.LBAtoms,i)==true{
			newbasis=fmt.Sprintf("newgto \"%s\" end",Q.LowBasis)
			}
		fmt.Fprintf(file,"%-2s  %8.3f%8.3f%8.3f %s\n",atoms.Atom(i).Symbol, coords.Get(i,0), coords.Get(i,2), coords.Get(i,2),newbasis)	
		}
	fmt.Fprintf(file,"*\n")
	return nil
	}
	
//Run runs the command given by the string O.command
//it waits or not for the result depending of 
func (O *OrcaRunner) Run(wait bool) (err error){
	command:=exec.Command("nohup",O.command,fmt.Sprintf("%s.inp",O.inputname))
	if wait==true{
		err=command.Run()
		}else{
		err=command.Start()	
		}
	return err
	}	
	
	
	
	
//buildCConstraints transforms the list of cartesian constrains in the QMCalc structre
//into a string with ORCA-formatted cartesian constraints
func (O *OrcaRunner) buildCConstraints(C []int) string{
	if C==nil{
		return "\n" //no constraints
		}
	constraints:=make([]string,len(C)+3)
	constraints[0]="%geom Constraints\n"
	for key,val:=range(C){
		constraints[key+1]=fmt.Sprintf("         {C %d C}\n",val)
		}
	last:=len(constraints)-1
	constraints[last-1]="         end\n"
	constraints[last]=" end\n"
	final:=strings.Join(constraints,"")
	return final
	}




var orcaSCFTight = map[int] string {
      1: "",
	  2: "TightSCF",
      3: "VeryTightSCF",
}

var orcaSCFConv = map[int] string {
      1: "",
	  2: "SlowConv",
      3: "VerySlowConv",
}


var orcaDisp = map[string] string {
	 "nodisp":"",
      "D":  "VDW04",
	  "D2": "VDW06",
      "D3": "VDW10",
}



