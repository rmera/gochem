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

package qm

import "os"
import "strings"
import "fmt"
import "runtime"
import "os/exec"
import "strconv"
import "github.com/rmera/gochem"

//Note that the default methods and basis vary with each program, and even
//for a given program they are NOT considered part of the API, so they can always change.
type OrcaHandle struct {
	defmethod   string
	defbasis    string
	defauxbasis string
	previousMO  string
	bsse        string
	command     string
	inputname   string
	nCPU        int
}

func NewOrcaHandle() *OrcaHandle {
	run := new(OrcaHandle)
	run.SetDefaults()
	return run
}

//OrcaHandle methods

//Sets the number of CPU to be used
func (O *OrcaHandle) SetnCPU(cpu int) {
	O.nCPU = cpu
}

func (O *OrcaHandle) SetName(name string) {
	O.inputname = name
}

func (O *OrcaHandle) SetCommand(name string) {
	O.command = name
}

func (O *OrcaHandle) SetMOName(name string) {
	O.previousMO = name
}

/*Sets defaults for ORCA calculation. Default is a single-point at
revPBE/def2-SVP with RI, and all the available CPU with a max of
8. The ORCA command is set to $ORCA_PATH/orca, at least in
unix.*/
func (O *OrcaHandle) SetDefaults() {
	O.defmethod = "revPBE"
	O.defbasis = "def2-SVP"
	O.defauxbasis = "def2-SVP/J"
	O.command = os.ExpandEnv("${ORCA_PATH}/orca")
	if O.command == "/orca" { //if ORCA_PATH was not defined
		O.command = "./orca"
	}
	cpu := runtime.NumCPU()
	O.nCPU = cpu

}

//BuildInput builds an input for ORCA based int the data in atoms, coords and C.
//returns only error.
func (O *OrcaHandle) BuildInput(coords *chem.VecMatrix, atoms chem.ReadRef, Q *Calc) error {
	//Only error so far
	if atoms == nil || coords == nil {
		return fmt.Errorf("Missing charges or coordinates")
	}
	if Q.Basis == "" {
		fmt.Fprintf(os.Stderr, "no basis set assigned for ORCA calculation, will used the default %s, \n", O.defbasis)
		Q.Basis = O.defbasis
	}
	if Q.Method == "" {
		fmt.Fprintf(os.Stderr, "no method assigned for ORCA calculation, will used the default %s, \n", O.defmethod)
		Q.Method = O.defmethod
		Q.auxColBasis = "" //makes no sense for pure functional
		Q.auxBasis = fmt.Sprintf("%s/J", Q.Basis)
	}

	//Set RI or RIJCOSX if needed
	ri := ""
	if Q.RI && Q.RIJ {
		return fmt.Errorf("RI and RIJ cannot be activated at the same time")
	}
	if Q.RI {
		Q.auxBasis = Q.Basis + "/J"
		//	if !strings.Contains(Q.Others," RI "){
		ri = "RI"
	}
	if Q.RIJ {
		Q.auxBasis = Q.Basis + "/J"
		Q.auxColBasis = Q.Basis + "/C"
		//	if !strings.Contains(Q.Others,"RIJCOSX"){
		ri = "RIJCOSX"
	}

	disp := "D3"
	if Q.Disperssion != "" {
		disp = orcaDisp[Q.Disperssion]
	}
	opt := ""
	if Q.Optimize == true {
		opt = "Opt"
	}
	//If this flag is set we'll look for a suitable MO file.
	//If not found, we'll just use the default ORCA guess
	hfuhf := "RHF"
	if atoms.Multi() != 1 {
		hfuhf = "UHF"
	}
	moinp := ""
	if Q.OldMO == true {
		dir, _ := os.Open("./")     //This should always work, hence ignoring the error
		files, _ := dir.Readdir(-1) //Get all the files.
		for _, val := range files {
			if O.previousMO != "" {
				break
			}
			if val.IsDir() == true {
				continue
			}
			name := val.Name()
			if strings.Contains(".gbw", name) {
				O.previousMO = name
				break
			}
		}
		if O.previousMO != "" {
			Q.Guess = "MORead"
			moinp = fmt.Sprintf("%%moinp \"%s\"\n", O.previousMO)
		} else {
			moinp = ""
			Q.Guess = "" //The default guess
		}
	}
	tight := "TightSCF"
	if Q.SCFTightness != 0 {
		tight = orcaSCFTight[Q.SCFTightness]
	}
	conv := ""
	if Q.SCFConvHelp == 0 {
		//The default for this is nothing for RHF and SlowConv for UHF
		if atoms.Multi() > 1 {
			conv = "SlowConv"
		}
	} else {
		conv = orcaSCFConv[Q.SCFConvHelp]
	}
	pal := ""
	palbig := ""
	if O.nCPU > 1 {
		if O.nCPU > 8 {
			palbig = fmt.Sprintf("%%pal nprocs %d\n   end\n", O.nCPU) //fmt.Fprintf(os.Stderr, "CPU number of %d for ORCA calculations currently not supported, maximun 8", O.nCPU)
		} else {
			pal = fmt.Sprintf("PAL%d", O.nCPU)
		}
	}
	grid := ""
	if Q.Grid > 0 && Q.Grid <= 9 {
		final := ""
		if Q.Grid > 3 {
			final = "NoFinalGrid"
		}
		grid = fmt.Sprintf("Grid%d %s", Q.Grid, final)
	}
	var err error
	var bsse string
	if bsse, err = O.buildgCP(Q); err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
	}
	if Q.Method == "HF-3c" { //This method includes its own basis sets and corrections, so previous choices are overwritten.
		Q.Basis = ""
		Q.auxBasis = ""
		Q.auxColBasis = ""
		Q.Guess = ""
		bsse = ""
		disp = ""
	}
	MainOptions := []string{"!", hfuhf, Q.Method, Q.Basis, Q.auxBasis, Q.auxColBasis, tight, disp, conv, Q.Guess, opt, Q.Others, grid, ri, bsse, pal, "\n"}
	mainline := strings.Join(MainOptions, " ")
	constraints := O.buildCConstraints(Q.CConstraints)
	iconstraints, err := O.buildIConstraints(Q.IConstraints)
	if err != nil {
		return err
	}
	cosmo := ""
	if Q.Dielectric > 0 {
		cosmo = fmt.Sprintf("%%cosmo epsilon %1.0f\n        refrac 1.30\n        end\n", Q.Dielectric)
	}
	mem := ""
	if Q.Memory != 0 {
		mem = fmt.Sprintf("%%MaxCore %d\n", Q.Memory)
	}
	ElementBasis := ""
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
	//Now lets write the thing
	if O.inputname == "" {
		O.inputname = "gochem"
	}
	file, err := os.Create(fmt.Sprintf("%s.inp", O.inputname))
	if err != nil {
		return err
	}
	defer file.Close()
	_, err = fmt.Fprint(file, mainline)
	//With this check its assumed that the file is ok.
	if err != nil {
		return err
	}
	fmt.Fprint(file, palbig)
	fmt.Fprint(file, moinp)
	fmt.Fprint(file, mem)
	fmt.Fprint(file, constraints)
	fmt.Fprint(file, iconstraints)
	fmt.Fprint(file, ElementBasis)
	fmt.Fprint(file, cosmo)
	fmt.Fprint(file, "\n")
	//Now the type of coords, charge and multiplicity
	fmt.Fprintf(file, "* xyz %d %d\n", atoms.Charge(), atoms.Multi())
	//now the coordinates
	//	fmt.Println(atoms.Len(), coords.Rows()) ///////////////
	for i := 0; i < atoms.Len(); i++ {
		newbasis := ""
		if isInInt(Q.HBAtoms, i) == true {
			newbasis = fmt.Sprintf("newgto \"%s\" end", Q.HighBasis)
		} else if isInInt(Q.LBAtoms, i) == true {
			newbasis = fmt.Sprintf("newgto \"%s\" end", Q.LowBasis)
		}
		//	fmt.Println(atoms.Atom(i).Symbol)
		fmt.Fprintf(file, "%-2s  %8.3f%8.3f%8.3f %s\n", atoms.Atom(i).Symbol, coords.At(i, 0), coords.At(i, 1), coords.At(i, 2), newbasis)
	}
	fmt.Fprintf(file, "*\n")
	return nil
}

//Run runs the command given by the string O.command
//it waits or not for the result depending on wait.
//Not waiting for results works
//only for unix-compatible systems, as it uses bash and nohup.
func (O *OrcaHandle) Run(wait bool) (err error) {

	if wait == true {
		out, err := os.Create(fmt.Sprintf("%s.out", O.inputname))
		if err != nil {
			return err
		}
		defer out.Close()
		command := exec.Command(O.command, fmt.Sprintf("%s.inp", O.inputname))
		command.Stdout = out
		err = command.Run()

	} else {
		command := exec.Command("sh", "-c", "nohup "+O.command+fmt.Sprintf(" %s.inp > %s.out &", O.inputname, O.inputname))
		err = command.Start()
	}
	return err
}

//buildIConstraints transforms the list of cartesian constrains in the QMCalc structre
//into a string with ORCA-formatted internal constraints.
func (O *OrcaHandle) buildIConstraints(C []*IConstraint) (string, error) {
	if C == nil {
		return "\n", nil //no constraints
	}
	constraints := make([]string, len(C)+3)
	constraints[0] = "%geom Constraints\n"
	for key, val := range C {

		if iConstraintOrder[val.Class] != len(val.CAtoms) {
			return "", fmt.Errorf("Internal constraint ill-formated")
		}

		var temp string
		if val.Class == 'B' {
			temp = fmt.Sprintf("         {B %d %d %2.3f C}\n", val.CAtoms[0], val.CAtoms[1], val.Val)
		} else if val.Class == 'A' {
			temp = fmt.Sprintf("         {A %d %d %d %2.3f C}\n", val.CAtoms[0], val.CAtoms[1], val.CAtoms[2], val.Val)
		} else if val.Class == 'D' {
			temp = fmt.Sprintf("         {D %d %d %d %d %2.3f C}\n", val.CAtoms[0], val.CAtoms[1], val.CAtoms[2], val.CAtoms[3], val.Val)
		}
		constraints[key+1] = temp
	}
	last := len(constraints) - 1
	constraints[last-1] = "         end\n"
	constraints[last] = " end\n"
	final := strings.Join(constraints, "")
	return final, nil
}

var iConstraintOrder = map[byte]int{
	'B': 2,
	'A': 3,
	'D': 4,
}

//buildCConstraints transforms the list of cartesian constrains in the QMCalc structre
//into a string with ORCA-formatted cartesian constraints
func (O *OrcaHandle) buildCConstraints(C []int) string {
	if C == nil {
		return "\n" //no constraints
	}
	constraints := make([]string, len(C)+3)
	constraints[0] = "%geom Constraints\n"
	for key, val := range C {
		constraints[key+1] = fmt.Sprintf("         {C %d C}\n", val)
	}
	last := len(constraints) - 1
	constraints[last-1] = "         end\n"
	constraints[last] = " end\n"
	final := strings.Join(constraints, "")
	return final
}

//Only DFT is supported. Also, only Karlsruhe's basis sets. If you are using Pople's,
//let us know so we can send a mission to rescue you from the sixties :-)
func (O *OrcaHandle) buildgCP(Q *Calc) (string, error) {
	ret := ""
	var err error
	if strings.ToLower(Q.BSSE) == "gcp" {
		switch strings.ToLower(Q.Basis) {
		case "def2-svp":
			ret = "GCP(DFT/SVP)"
		case "def2-tzvp":
			ret = "GCP(DFT/TZ)"
		case "def2-sv(p)":
			ret = "GCP(DFT/SV(P))"
		default:
			err = fmt.Errorf("Method/basis combination for gCP unavailable, will skip the correction")
		}
	}
	return ret, err
}

var orcaSCFTight = map[int]string{
	0: "",
	1: "TightSCF",
	2: "VeryTightSCF",
}

var orcaSCFConv = map[int]string{
	0: "",
	1: "SlowConv",
	2: "VerySlowConv",
}

var orcaDisp = map[string]string{
	"nodisp": "",
	"D3OLD":  "VDW10", //compatibility with ORCA 2.9
	"D2":     "D2",
	"D3BJ":   "D3BJ",
	"D3bj":   "D3BJ",
	"D3":     "D3ZERO",
	"D3ZERO": "D3ZERO",
	"D3Zero": "D3ZERO",
	"D3zero": "D3ZERO",
	"VV10":   "NL", //for these methods only the default integration grid is supported.
	"SCVV10": "SCNL",
	"NL":     "NL",
	"SCNL":   "SCNL",
}

/*Reads the latest geometry from an ORCA optimization. Returns the
  geometry or error. Returns the geometry AND error if the geometry read
  is not the product of a correctly ended ORCA calculation. In this case
  the error is "probable problem in calculation"*/
func (O *OrcaHandle) OptimizedGeometry(atoms chem.Ref) (*chem.VecMatrix, error) {
	var err error
	geofile := fmt.Sprintf("%s.xyz", O.inputname)
	//Here any error of orcaNormal... or false means the same, so the error can be ignored.
	if trust := O.orcaNormalTermination(); !trust {
		err = fmt.Errorf("Probable problem in calculation")
	}
	//This might not be super efficient but oh well.
	mol, err1 := chem.XYZRead(geofile)
	if err1 != nil {
		return nil, err1
	}
	return mol.Coords[0], err //returns the coords, the error indicates whether the structure is trusty (normal calculation) or not
}

//Gets the energy of a previous Orca calculations.
//Returns error if problem, and also if the energy returned that is product of an
//abnormally-terminated ORCA calculation. (in this case error is "Probable problem
//in calculation")
func (O *OrcaHandle) Energy() (float64, error) {
	err := fmt.Errorf("Probable problem in calculation")
	f, err1 := os.Open(fmt.Sprintf("%s.out", O.inputname))
	if err1 != nil {
		return 0, err1
	}
	defer f.Close()
	f.Seek(-1, 2) //We start at the end of the file
	energy := 0.0
	var found bool
	for i := 0; ; i++ {
		line, err1 := getTailLine(f)
		if err1 != nil {
			return 0.0, err1
		}
		if strings.Contains(line, "**ORCA TERMINATED NORMALLY**") {
			err = nil
		}
		if strings.Contains(line, "FINAL SINGLE POINT ENERGY") {
			splitted := strings.Fields(line)
			energy, err1 = strconv.ParseFloat(splitted[4], 64)
			if err1 != nil {
				return 0.0, err1
			}
			found = true
			break
		}
	}
	if !found {
		return 0.0, fmt.Errorf("Output does not contain energy")
	}
	return energy * chem.H2Kcal, err
}

//Gets previous line of the file f
func getTailLine(f *os.File) (line string, err error) {
	var i int64 = 1
	buf := make([]byte, 1)
	var ini int64 = -1
	for ; ; i++ {
		//move the pointer back one byte per cycle
		if _, err := f.Seek(-2, 1); err != nil {
			return "", err
		}
		if _, err := f.Read(buf); err != nil {
			return "", err
		}
		if buf[0] == byte('\n') && ini == -1 {
			ini = i
			break
		}
	}
	bufF := make([]byte, ini)
	f.Read(bufF)

	if _, err := f.Seek(int64(-1*(len(bufF))), 1); err != nil { //making up for the read
		return "", err
	}
	return string(bufF), nil
}

//This checks that an ORCA calculation has terminated normally
//I know this duplicates code, I wrote this one first and then the other one.
func (O *OrcaHandle) orcaNormalTermination() bool {
	var ini int64 = 0
	var end int64 = 0
	var first bool
	buf := make([]byte, 1)
	f, err := os.Open(fmt.Sprintf("%s.out", O.inputname))
	if err != nil {
		return false
	}
	defer f.Close()
	var i int64 = 1
	for ; ; i++ {
		if _, err := f.Seek(-1*i, 2); err != nil {
			return false
		}
		if _, err := f.Read(buf); err != nil {
			return false
		}
		if buf[0] == byte('\n') && first == false {
			first = true
		} else if buf[0] == byte('\n') && end == 0 {
			end = i
		} else if buf[0] == byte('\n') && ini == 0 {
			ini = i
			break
		}

	}
	f.Seek(-1*ini, 2)
	bufF := make([]byte, ini-end)
	f.Read(bufF)
	if strings.Contains(string(bufF), "**ORCA TERMINATED NORMALLY**") {
		return true
	}
	return false
}
