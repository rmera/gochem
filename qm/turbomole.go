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

package qm

import "os"
import "io"
import "strings"
import "log"

import "strconv"
import "bufio"
import "fmt"
import "os/exec"
import "github.com/rmera/gochem"

//This imlpementation supports only singlets and doublets.
type TMHandle struct {
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
func NewTMHandle() *TMHandle {
	run := new(TMHandle)
	run.SetDefaults()
	return run
}

//TMHandle methods

//Just to satisfy the interface. It does nothing
func (O *TMHandle) SetnCPU(cpu int) {
	//It does nothing! :-D
}

//This set the name of the subdirectory, in the current directory
//where the calculation will be ran
func (O *TMHandle) SetName(name string) {
	O.inputname = name

}

//In TM the command is set according to the method. I just assume a normal installation.
//This method doesnt do anything.
func (O *TMHandle) SetCommand(name string) {
	//Does nothing again
}

//Sets some defaults for TMHandle. default is an optimization at
//  TPSS-D3 / def2-SVP
func (O *TMHandle) SetDefaults() {
	O.defmethod = "tpss"
	O.defbasis = "def2-SVP"
	O.defauxbasis = "def2-SVP"
	O.command = "ridft"
	O.inputname = "gochemturbo"

}

//Adds all the strings in toapend to the control file, just before the $symmetry keyword
func (O *TMHandle) addToControl(toappend []string, Q *Calc) error {
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

func (O *TMHandle) addCosmo(epsilon float64) error {
	//The ammount of newlines is wrong, must fix
	cosmostring := "" //a few newlines before the epsilon
	if epsilon == 0 {
		return nil
	}
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

func (O *TMHandle) addBasis(basisOrEcp string, basiselems []string, basis, defstring string) string {
	if basiselems == nil { //no atoms to add basis to, do nothing
		return defstring
	}
	for _, elem := range basiselems {
		defstring = fmt.Sprintf("%s%s \"%s\" %s\n", defstring, basisOrEcp, strings.ToLower(elem), basis)
	}
	return defstring
}

//modifies the coord file such as to freeze the atoms in the slice frozen.
func (O *TMHandle) addFrozen(frozen []int) error {
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
//Note that at this point the interface does not support multiplicities different from 1 and 2.
//The number in atoms is simply ignored.
func (O *TMHandle) BuildInput(coords *chem.VecMatrix, atoms chem.ReadRef, Q *Calc) error {
	err := os.Mkdir(O.inputname, os.FileMode(0755))
	for i := 0; err != nil; i++ {
		if strings.Contains(err.Error(), "file exists") {
			O.inputname = fmt.Sprintf("%s%d", O.inputname, i)
			err = os.Mkdir(O.inputname, os.FileMode(0755))
		} else {
			return err
		}
	}
	_ = os.Chdir(O.inputname)
	defer os.Chdir("..")
	//Set the coordinates in a slightly stupid way.
	chem.XYZWrite("file.xyz", coords, atoms)
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
	defstring = O.addBasis("b", Q.LBElements, Q.LowBasis, defstring)
	defstring = O.addBasis("ecp", Q.ECPElements, Q.ECP, defstring)
	defstring = O.addBasis("b", Q.ECPElements, Q.ECP, defstring)      //we set a basis set compatible with the ECP. In TM they share the same name
	defstring = O.addBasis("b", Q.HBElements, Q.HighBasis, defstring) //The high basis will override the ECP basis, which can be rather small. Use under your own risk.
	defstring = defstring + "\n\n\n\n*\n"
	//This is because the %$$#%^ define interface ask some $#%&&# questions in the eht setup when encounters some atoms.
	//so i have to add an additional newline for each of these types. So far I know only that copper causes this.
	stupid := ""
	stupidatoms := "Cu" //if you want to add more stupid atoms jsut add then to the string: "Cu Zn"
	for i := 0; i < atoms.Len(); i++ {
		if stupidatoms == "" {
			break
		}
		if strings.Contains(stupidatoms, atoms.Atom(i).Symbol) {
			stupidatoms = strings.Replace(stupidatoms, atoms.Atom(i).Symbol, "", -1)
			stupid = stupid + "\n"
		}
	}
	defstring = fmt.Sprintf("%seht\n%s\n%d\n\n", defstring, stupid, atoms.Charge())
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
	log.Println(defstring)
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
	err = O.addCosmo(Q.Dielectric)
	return err
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
	"D":      "$olddisp",
	"D2":     "$disp2",
	"D3":     "$disp3",
}

//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//This is a Unix-only function.
func (O *TMHandle) Run(wait bool) (err error) {
	os.Chdir(O.inputname)
	defer os.Chdir("..")
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

//Energy returns the energy from the corresponding calculation, in kcal/mol.
func (O *TMHandle) Energy() (float64, error) {
	os.Chdir(O.inputname)
	defer os.Chdir("..")
	f, err := os.Open("energy")
	if err != nil {
		return 0, err
	}
	defer f.Close()
	fio := bufio.NewReader(f)
	line, err := getSecondToLastLine(fio)
	if err != nil {
		return 0, err
	}
	en := strings.Fields(line)[1]
	energy, err := strconv.ParseFloat(en, 64)
	return energy * chem.H2Kcal, err
}

//OptimizedGeometry returns the coordinates for the optimized structure.
func (O *TMHandle) OptimizedGeometry(atoms chem.Ref) (*chem.VecMatrix, error) {
	os.Chdir(O.inputname)
	defer os.Chdir("..")
	x2t := exec.Command("t2x")
	stdout, err := x2t.StdoutPipe()
	if err != nil {
		return nil, fmt.Errorf("Unable to run t2x: %s", err.Error())
	}
	if err := x2t.Start(); err != nil {
		return nil, fmt.Errorf("Unable to run t2x: %s", err.Error())
	}
	xyz := bufio.NewReader(stdout)
	mol, err := chem.XYZBufIORead(xyz)
	if err != nil {
		return nil, err
	}
	return mol.Coords[len(mol.Coords)-1], nil

}

//Gets the second to last line in a turbomole energy file given as a bufio.Reader.
//expensive on the CPU but rather easy on the memory, as the file is read line by line.
func getSecondToLastLine(f *bufio.Reader) (string, error) {
	prevline := ""
	line := ""
	var err error
	for {
		line, err = f.ReadString('\n')
		if err != nil {
			break
		}
		if !strings.Contains(line, "$end") {
			prevline = line
		} else {
			break
		}
	}
	return prevline, err
}
