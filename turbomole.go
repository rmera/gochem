/*
 * turbomole.go, part of gochem.
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

//import "strconv"
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

//Sets some defaults for TMRunner. default is an optimization at
//  TPSS-D3 / def2-SVP
func (O *TMRunner) SetDefaults() {
	O.defmethod = "tpss"
	O.defbasis = "def2-SVP"
	O.defauxbasis = "def2-SVP"
	O.command = "jobex -c 100 -ri > jobex.out"

}

//Adds all the strings in toapend to the control file, just before the $symmetry keyword
func (O *TMRunner) addToControl(toappend []string, Q *QMCalc) error {
	f, err := os.Open("control")
	if err != nil {
		return err
	}
	lines := make([]string, 0, 200) //200 is just a guess for the number of lines in the control file
	c := bufio.NewReader(f)
	for err == nil {
		var line string
		line, err = c.ReadString('\n')
		lines = append(lines, line)
	}
	f.Close() //I cant defer it because I need it closed now.
	out, err := os.Create("control")
	if err != nil {
		return err
	}
	defer out.Close()
	var k string
	for _, i := range lines {
		k = i //May not be too efficient
		if strings.Contains(i, "$symmetry") {
			for _, j := range toappend {
				if _, err := fmt.Fprintf(out, j+"\n"); err != nil {
					return err
				}
			}
		}
		if Q.SCFConvHelp >= 1 {
			if strings.Contains(k, "$scfiterlimit") {
				k = "$scfiterlimit   100\n"
			}
			if strings.Contains(k, "$scfdamp") {
				k = "$scfdamp start=10 step=0.005 min=0.5\n"
			}
		}
		if _, err := fmt.Fprintf(out, k); err != nil {
			return err
		}
	}
	return nil
}

func (O *TMRunner) addCosmo(epsilon float64) error {
	//The ammount of newlines is wrong, must fix
	cosmostring := "" //a few newlines before the epsilon
	cosmostring = fmt.Sprintf("%s%3.1f\n\n\n\n\n\n\n\nr all b\n*\n\n", cosmostring, epsilon)
	def := exec.Command("cosmoprep")
	pipe, err := def.StdinPipe()
	if err != nil {
		return fmt.Errorf("Unable to run cosmoprep: %s", err.Error())
	}
	defer pipe.Close()
	pipe.Write([]byte(cosmostring))
	if err := def.Run(); err != nil {
		return fmt.Errorf("Unable to run run cosmoprep: %s", err.Error())
	}
	return nil

}

func (O *TMRunner) addBasis(basiselems []string, basis, defstring string) string {
	if basiselems == nil { //no atoms to add basis to, do nothing
		return defstring
	}
	for _, elem := range basiselems {
		defstring = fmt.Sprintf("%sb \"%s\" %s\n", defstring, strings.ToLower(elem), basis)
	}
	return defstring
}

//modifies the coord file such as to freeze the atoms in the slice frozen.
func (O *TMRunner) addFrozen(frozen []int) error {
	f, err := os.Open("coord")
	if err != nil {
		return err
	}
	lines := make([]string, 0, 200) //200 is just a guess for the number of lines in the coord file
	c := bufio.NewReader(f)
	for err == nil {
		var line string
		line, err = c.ReadString('\n')
		lines = append(lines, line)
	}
	f.Close() //I cant defer it because I need it closed now.
	out, err := os.Create("coord")
	if err != nil {
		return err
	}
	defer out.Close()
	for key, i := range lines {
		if isInInt(frozen, key-1) {
			j := strings.Replace(i, "\n", " f\n", -1)
			if _, err := fmt.Fprintf(out, j); err != nil {
				return err
			}
		} else {
			if _, err := fmt.Fprintf(out, i); err != nil {
				return err
			}
		}
	}
	return nil
}

func copy2pipe(pipe io.ReadCloser, file *os.File, end chan bool) {
	io.Copy(file, pipe)
	end <- true
}

//BuildInput builds an input for TM based int the data in atoms, coords and C.
//returns only error.
func (O *TMRunner) BuildInput(atoms Ref, coords *CoordMatrix, Q *QMCalc) error {
	//Set the coordinates in a slightly stupid way.
	XyzWrite("file.xyz", atoms, coords)
	x2t := exec.Command("x2t", "file.xyz")
	stdout, err := x2t.StdoutPipe()
	if err != nil {
		return fmt.Errorf("Unable to run x2t: %s", err.Error())
	}
	coord, err := os.Create("coord")
	if err != nil {
		return fmt.Errorf("Unable to run x2t: %s", err.Error())
	}
	if err := x2t.Start(); err != nil {
		return fmt.Errorf("Unable to run x2t: %s", err.Error())
	}
	//	var end chan bool
	//	go copy2pipe(stdout, coord, end)
	//	<-end
	io.Copy(coord, stdout)
	coord.Close() //not defearable
	defstring := "\n\na coord\n*\nno\n"
	if atoms == nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
	}
	if Q.Basis == "" {
		fmt.Fprintf(os.Stderr, "no basis set assigned for TM calculation, will used the default %s, \n", O.defbasis)
		Q.Basis = O.defbasis
	}
	defstring = defstring + "b all " + Q.Basis + "\n"
	defstring = O.addBasis(Q.HBElements, Q.HighBasis, defstring)
	defstring = O.addBasis(Q.LBElements, Q.LowBasis, defstring)
	defstring = defstring + "\n\n\n\n*\n"
	defstring = fmt.Sprintf("%seht\n\n\n%d\n\n", defstring, atoms.Charge())
	method, ok := tMMethods[Q.Method]
	if !ok {
		fmt.Fprintf(os.Stderr, "no method assigned for TM calculation, will used the default %s, \n", O.defmethod)
		Q.Method = O.defmethod
		Q.RI = true
	} else {
		Q.Method = method
	}
	//We only support HF and DFT
	//O.command = "dscf"
	if Q.Method != "hf" {
		grid := ""
		if Q.Grid != 0 && Q.Grid <= 7 {
			grid = fmt.Sprintf("grid\n m%d\n", Q.Grid)
		}
		defstring = defstring + "dft\non\nfunc " + Q.Method + "\n" + grid + "*\n"
		if Q.RI {
			mem := 500
			if Q.Memory != 0 {
				mem = Q.Memory
			}
			defstring = fmt.Sprintf("%sri\non\nm %d\n*\n", defstring, mem)
			O.command = "ridft"
		}
	}
	defstring = defstring + "*\n"
	fmt.Println(defstring) ////////////////
	def := exec.Command("define")
	pipe, err := def.StdinPipe()
	if err != nil {
		return fmt.Errorf("Unable to run define: %s", err.Error())
	}
	defer pipe.Close()
	pipe.Write([]byte(defstring))
	if err := def.Run(); err != nil {
		return fmt.Errorf("Unable to run define: %s", err.Error())
	}

	if Q.Optimize {
		O.command = "jobex"
		if Q.RI {
			O.command = O.command + " -c 200 -ri"
		} else {
			O.command = O.command + " -c 200"
		}
	}
	//Now modify control
	args := make([]string, 1, 2)
	args[0], ok = tMDisp[Q.Disperssion]
	if !ok {
		fmt.Fprintf(os.Stderr, "Dispersion correction requested not supported, will used the default: D3, \n")
		args[0] = "$disp3"
	}
	if Q.Gimic {
		O.command = "mpshift"
		args = append(args, "$gimic")
	}
	if err := O.addToControl(args, Q); err != nil {
		return err
	}
	//set the frozen atoms (only cartesian constraints are supported)
	if err := O.addFrozen(Q.CConstraints); err != nil {
		return err
	}
	//Finally the cosmo business.
	return O.addCosmo(Q.Dielectric)
}

var tMMethods = map[string]string{
	"HF":     "hf",
	"hf":     "hf",
	"b3lyp":  "b3-lyp",
	"B3LYP":  "b3-lyp",
	"b3-lyp": "b3-lyp",
	"PBE":    "pbe",
	"pbe":    "pbe",
	"TPSS":   "tpss",
	"TPSSh":  "tpssh",
	"tpss":   "tpss",
	"tpssh":  "tpssh",
	"BP86":   "b-p",
	"b-p":    "b-p",
}

var tMDisp = map[string]string{
	"":       "",
	"nodisp": "",
	"D":      "$disp",
	"D2":     "$disp2",
	"D3":     "$disp3",
}

//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//This is a Unix-only function.
func (O *TMRunner) Run(wait bool) (err error) {
	filename := strings.Fields(O.command)
	fmt.Println("nohup " + O.command + " > " + filename[0] + ".out")
	command := exec.Command("sh", "-c", "nohup "+O.command+" >"+filename[0]+".out")
	if wait == true {
		err = command.Run()
	} else {
		err = command.Start()
	}
	return err
}

//GetEnergy is NOT working yet
func (O *TMRunner) GetEnergy() (float64, error) {
	return 0, nil
}

//GetGeometry is NOT working yet
func (O *TMRunner) GetGeometry(atoms Ref) (*matrix.DenseMatrix, error) {
	return nil, nil

}
