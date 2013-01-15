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
import "io"
import "strings"
import "strconv"
import "bufio"
import "fmt"
import "github.com/skelterjohn/go.matrix"
import "os/exec"


type TMRunner struct {
	defmethod   string
	defbasis    string
	defauxbasis string
	previousMO  string
	command     string
	inputname   string
	gimic       bool
}

//Creates and initialized a new instance of TMRuner, with values set
//to its defaults.
func MakeTMRunner() *TMRunner {
	run := new(TMRunner)
	run.SetDefaults()
	return run
}

//TMRunner methods


//Just to satisfy the interface. It does nothing
func (O *TMRunner) SetnCPU(cpu int) {
	//It does nothing! :-D
}

//This method does nothing. Just there to satisfy the interface
//In TM you dont quite get to choose names. We wont here at least
func (O *TMRunner) SetName(name string) {
}

//In TM the command is set according to the method. I just assume a normal installation.
//This method doesnt do anything.
func (O *TMRunner) SetCommand(name string) {
	//Does nothing again
}

/*Sets some defaults for TMRunner. default is an optimization at
  PM6-DH2X It tries to locate TM2009 according to the
  $TM_LICENSE environment variable, which might only work in UNIX. 
  If other system or using TM2012 the command Must be set with the
  SetCommand function. */
func (O *TMRunner) SetDefaults() {
	O.defmethod = "tpss"
	O.defbasis = "def2-SVP" 
	O.defauxbasis = "def2-SVP"
	O.command = "ridft"
	
}


func (O *TMRunner) addToControl(toappend []string) error{
	fmt.Println("running add2control")
	f, err:=os.OpenFile("control", os.O_APPEND, 0666)
	if err!=nil{
		return err
	}
	defer f.Close()
	fmt.Println("opened!")
	for _,val:=range(toappend){
		fmt.Println("writing!")
		if _, err := io.WriteString(f, val+"\n");err!=nil{
			return err
		} 
		fmt.Println("wrote!")	
	}
	return nil
}

//BuildInput builds an input for TM based int the data in atoms, coords and C.
//returns only error.
func (O *TMRunner) BuildInput(atoms Ref, coords *matrix.DenseMatrix, Q *QMCalc) error {	
	//Set the coordinates in a slightly stupid way.
	XyzWrite(atoms,coords,"file.xyz")
	x2t:=exec.Command("x2t","file.xyz")
	stdout, err := x2t.StdoutPipe()
	if err != nil {
		return fmt.Errorf("Unable to run x2t: %s",err.Error())
	}
    coord,err:=os.Create("coord")
    if err!=nil{
		return fmt.Errorf("Unable to run x2t: %s",err.Error())
		}
	defer coord.Close()
    go io.Copy(coord, stdout) 
	if err:=x2t.Run();err!=nil{
		return fmt.Errorf("Unable to run x2t: %s",err.Error())
	}	
	
	fmt.Println("READY!")
	defstring:="\n\na coord\n*\nno\n"
	if atoms == nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
	}
	if Q.Basis == "" {
		fmt.Fprintf(os.Stderr, "no basis set assigned for TM calculation, will used the default %s, \n", O.defbasis)
		defstring=defstring+"b a "+Q.Basis+"\n*\n"
	}
	defstring=defstring+fmt.Sprintf("eht\n\n\n%d\n\n", atoms.Charge())
	if method,ok:=tMMethods[Q.Method]; !ok{
		fmt.Fprintf(os.Stderr, "no method assigned for TM calculation, will used the default %s, \n", O.defmethod)
		Q.Method = O.defmethod
		Q.RI = true
	}else{
	Q.Method=method	
	}
	//We only support HF and DFT
	O.command="dscf"
	if Q.Method!="hf" {
		defstring=defstring+"dft\non\nfunc "+Q.Method+"\n*\n"	
		if Q.RI{
			defstring=defstring+"ri\non\nm 500\n*\n"
			O.command="ridft"
			}
	}
	defstring=defstring+"*\n"
	
	def:=exec.Command("define")
	pipe,err:=def.StdinPipe()
	if err!=nil{
		return fmt.Errorf("Unable to run define: %s",err.Error())
		}
	defer pipe.Close()
	pipe.Write([]byte(defstring))
	if err:=def.Run();err!=nil{
		return fmt.Errorf("Unable to run define: %s",err.Error())
		}
		
	if Q.Optimize{
		O.command="jobex"
		if Q.RI{
			O.command=O.command+" -c 200 -ri"
		}else{
			O.command=O.command+" -c 200"
		}
	}
	//Now modify control
	args:=make([]string,1,2)
	args[0] = "$disp3"
	if Q.Disperssion != "" {
		args[0] = tMDisp[Q.Disperssion]	
	}
	if Q.Gimic{
		O.command="mpshift"
		args=append(args,"$gimic")
	}
	if err:=O.addToControl(args);err!=nil{
		return err
	}

/*	ElementBasis := ""
	if Q.HBElements != nil || Q.LBElements != nil {
		elementbasis := make([]string, 0, len(Q.HBElements)+len(Q.LBElements)+2)
		elementbasis = append(elementbasis, " %basis \n")
		for _, val := range Q.HBElements {
			elementbasis = append(elementbasis, fmt.Sprintf("  newgto %s \"%s\" end\n", val, Q.HighBasis))
		}
		for _, val := range Q.LBElements {
			elementbasis = append(elementbasis, fmt.Sprintf("  newgto %s \"%s\" end\n", val, Q.LowBasis))
		}
		elementbasis = append(elementbasis, "         end\n")
		ElementBasis = strings.Join(elementbasis, "")
	}
*/
	return nil
}

var tMMethods = map[string]string{
	"HF": "hf",
	"hf": "hf",
	"b3lyp": "b3-lyp",
	"B3LYP": "b3-lyp",
	"b3-lyp": "b3-lyp",
	"PBE": "pbe",
	"pbe": "pbe",
	"TPSS": "tpss",
	"TPSSh": "tpssh",
	"tpss": "tpss",
	"tpssh": "tpssh",
	"BP86": "b-p",
	"b-p": "b-p",
}


var tMDisp = map[string]string{
	"nodisp": "",
	"D":      "$disp",
	"D2":     "$disp2",
	"D3":     "$disp3",
}

//Run runs the command given by the string O.command
//it waits or not for the result depending of 
func (O *TMRunner) Run(wait bool) (err error) {
	command := exec.Command("nohup", O.command, fmt.Sprintf("%s.mop", O.inputname))
	if wait == true {
		err = command.Run()
	} else {
		err = command.Start()
	}
	return err
}

/*GetEnergy gets the last energy for a TM2009/2012 calculation by
  parsing the TM output file. Return error if fail. Also returns
  Error ("Probable problem in calculation")
  if there is a energy but the calculation didnt end properly*/
func (O *TMRunner) GetEnergy() (float64, error) {
	var err error
	var energy float64
	file, err := os.Open(fmt.Sprintf("%s.out", O.inputname))
	if err != nil {
		return 0, err
	}
	defer file.Close()
	out := bufio.NewReader(file)
	err = fmt.Errorf("TM Energy not found in %s", O.inputname)
	trust_radius_warning := false
	for {
		var line string
		line, err = out.ReadString('\n')
		if err != nil {
			break
		}
		if strings.Contains(line, "TRUST RADIUS NOW LESS THAN 0.00010 OPTIMIZATION TERMINATING") {
			trust_radius_warning = true
			continue
		}
		if strings.Contains(line, "TOTAL ENERGY") {
			splitted := strings.Fields(line)
			if len(splitted) < 4 {
				err = fmt.Errorf("Error reading energy from TM output file!")
				break
			}
			energy, err = strconv.ParseFloat(splitted[3], 64)
			if err != nil {
				break
			}
			energy = energy * eV2Kcalmol
			err = nil
			break
		}
	}
	if err != nil {
		return 0, err
	}
	if trust_radius_warning {
		err = fmt.Errorf("Probable problem in calculation")
	}
	return energy, err
}

/*Get Geometry reads the optimized geometry from a TM2009/2012 output. 
  Return error if fail. Returns Error ("Probable problem in calculation")
  if there is a geometry but the calculation didnt end properly*/
func (O *TMRunner) GetGeometry(atoms Ref) (*matrix.DenseMatrix, error) {
	var err error
	natoms := atoms.Len()
	coords := make([]float64, natoms*3, natoms*3) //will be used for return
	file, err := os.Open(fmt.Sprintf("%s.out", O.inputname))
	if err != nil {
		return nil, err
	}
	defer file.Close()
	out := bufio.NewReader(file)
	err = fmt.Errorf("TM Energy not found in %s", O.inputname)
	//some variables that will be changed/increased during the next for loop
	final_point := false //to see if we got to the right part of the file
	reading := false     //start reading
	i := 0
	errsl := make([]error, 3, 3)
	trust_radius_warning := false
	for {
		var line string
		line, err = out.ReadString('\n')
		if err != nil {
			break
		}

		if (!reading) && strings.Contains(line, "TRUST RADIUS NOW LESS THAN 0.00010 OPTIMIZATION TERMINATING") {
			trust_radius_warning = true
			continue
		}

		if !reading && strings.Contains(line, "FINAL  POINT  AND  DERIVATIVES") {
			final_point = true
			continue
		}
		if strings.Contains(line, "(ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)") && final_point {
			_, err = out.ReadString('\n')
			if err != nil {
				break
			}
			reading = true
			continue
		}
		if reading {
			//So far we dont check that there are not too many atoms in the TM output.
			if i*3+2 > natoms*3 {
				err = nil
				break
			}
			coords[i*3], errsl[0] = strconv.ParseFloat(strings.TrimSpace(line[22:35]), 64)
			coords[i*3+1], errsl[1] = strconv.ParseFloat(strings.TrimSpace(line[38:51]), 64)
			coords[i*3+2], errsl[2] = strconv.ParseFloat(strings.TrimSpace(line[54:67]), 64)
			i++
			err = parseErrorSlice(errsl)
			if err != nil {
				break
			}
		}
	}
	if err != nil {
		return nil, err
	}
	mcoords := matrix.MakeDenseMatrix(coords, natoms, 3)
	if trust_radius_warning {
		return mcoords, fmt.Errorf("Probable problem in calculation")
	}
	return mcoords, nil
}

