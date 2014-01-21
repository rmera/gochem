/*
 * files.go, part of gochem.
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

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

//PDB_read family

//A map for assigning mass to elements.
//Note that just common "bio-elements" are present
var symbolMass = map[string]float64{
	"H":  1.0,
	"C":  12.01,
	"O":  16.00,
	"N":  14.01,
	"P":  30.97,
	"S":  32.06,
	"Se": 78.96,
	"K":  39.1,
	"Ca": 40.08,
	"Mg": 24.30,
	"Cl": 35.45,
	"Na": 22.99,
	"Cu": 63.55,
	"Zn": 65.38,
	"Co": 58.93,
	"Fe": 55.84,
	"Mn": 54.94,
	"Si": 28.08,
}

//A map between 3-letters name for aminoacidic residues to the corresponding 1-letter names.
var three2OneLetter = map[string]byte{
	"SER": 'S',
	"THR": 'T',
	"ASN": 'N',
	"GLN": 'Q',
	"SEC": 'U', //Selenocysteine!
	"CYS": 'C',
	"GLY": 'G',
	"PRO": 'P',
	"ALA": 'A',
	"VAL": 'V',
	"ILE": 'I',
	"LEU": 'L',
	"MET": 'M',
	"PHE": 'F',
	"TYR": 'Y',
	"TRP": 'W',
	"ARG": 'R',
	"HIS": 'H',
	"LYS": 'K',
	"ASP": 'D',
	"GLU": 'E',
}

//This tries to guess a chemical element symbol from a PDB atom name. Mostly based on AMBER names.
//It only deals with some common bio-elements.
func symbolFromName(name string) (string, error) {
	symbol := ""
	if len(name) == 1 {
		symbol = name //should work
	} else if len(name) == 4 || name[0] == 'H' { //I thiiink only Hs can have 4-char names in amber.
		symbol = "H"
		//it name has more than one character but less than four.
	} else if name[0] == 'C' {
		if name[0:2] == "CU" {
			symbol = "Cu"
		} else if name == "CO" {
			symbol = "Co"
		} else if name == "CL" {
			symbol = "Cl"
		} else {
			symbol = "C"
		}
	} else if name[0] == 'N' {
		if name == "NA" {
			symbol = "Na"
		} else {
			symbol = "N"
		}
	} else if name[0] == 'O' {
		symbol = "O"
	} else if name[0] == 'P' {
		symbol = "P"
	} else if name[0] == 'S' {
		if name == "SE" {
			symbol = "Se"
		} else {
			symbol = "S"
		}
	} else if name[0:2] == "ZN" {
		symbol = "Zn"
	}
	if symbol == "" {
		return symbol, fmt.Errorf("Couldn't guess symbol from PDB name")
	}
	return symbol, nil
}

//Parses a valid ATOM or HETATM line of a PDB file, returns an Atom
// object with the info except for the coordinates and b-factors, which  are returned
// separately as an array of 3 float64 and a float64, respectively
func read_full_pdb_line(line string, read_additional bool, contlines int) (*Atom, []float64, float64, error) {
	err := make([]error, 7, 7) //accumulate errors to check at the end of the readed line.
	coords := make([]float64, 3, 3)
	atom := new(Atom)
	atom.Het = strings.HasPrefix(line, "HETATM") //this is called twice in the worst case. should fix
	atom.Id, err[0] = strconv.Atoi(strings.TrimSpace(line[6:12]))
	atom.Name = strings.TrimSpace(line[12:16])
	//PDB says that pos. 17 is for other thing but I see that is
	//used for residue name in many cases*/
	atom.Molname = line[17:20]
	atom.Molname1 = three2OneLetter[atom.Molname]
	atom.Chain = string(line[21])
	atom.Molid, err[1] = strconv.Atoi(strings.TrimSpace(line[22:30]))
	//Here we shouldn't need TrimSpace, but I keep it just in case someone
	// doesn's use all the fields when writting a PDB*/
	coords[0], err[2] = strconv.ParseFloat(strings.TrimSpace(line[30:38]), 64)
	coords[1], err[3] = strconv.ParseFloat(strings.TrimSpace(line[38:46]), 64)
	coords[2], err[4] = strconv.ParseFloat(strings.TrimSpace(line[46:54]), 64)
	atom.Occupancy, err[5] = strconv.ParseFloat(strings.TrimSpace(line[54:60]), 64)
	//If I try not to declare this and just use :=, I get an "expected identifier" error
	var bfactor float64
	bfactor, err[6] = strconv.ParseFloat(strings.TrimSpace(line[60:66]), 64)
	//we try to read the additional only if indicated and if it is there
	// In this part we don't catch errors. If something is missing we
	// just ommit it
	if read_additional && len(line) >= 80 {
		atom.Symbol = strings.TrimSpace(line[76:78])
		atom.Symbol = strings.Title(atom.Symbol)
		atom.Charge = float64(line[78]) //strconv.ParseFloat(strings.TrimSpace(line[78:78]),64)
		if strings.Contains(line[79:79], "-") {
			atom.Charge = -1.0 * atom.Charge
		}
	}

	//This part tries to guess the symbol from the atom name, if it has not been read
	//No error checking here, just fills symbol with the empty string the function returns
	if len(atom.Symbol) == 0 {
		atom.Symbol, _ = symbolFromName(atom.Name)
	}

	for i := range err {
		if err[i] != nil {
			//Here I should add the line number to the returned error.
			return nil, nil, 0, err[i]
		}
	}
	if atom.Symbol != "" {
		atom.Mass = symbolMass[atom.Symbol] //Not error checking
	}
	return atom, coords, bfactor, nil
}

/*Parses a PDB line if only the coordinates and bfactors are to be read*/
func read_onlycoords_pdb_line(line string, contlines int) ([]float64, float64, error) {
	coords := make([]float64, 3, 3)
	err := make([]error, 4, 4)
	var bfactor float64 //I dont get why I must declare this instead of using :=
	//I get an "expected identifier" error if I do so.
	coords[0], err[0] = strconv.ParseFloat(strings.TrimSpace(line[30:38]), 64)
	coords[1], err[1] = strconv.ParseFloat(strings.TrimSpace(line[38:46]), 64)
	coords[2], err[2] = strconv.ParseFloat(strings.TrimSpace(line[46:54]), 64)
	bfactor, err[3] = strconv.ParseFloat(strings.TrimSpace(line[60:66]), 64)
	//this will take care of any error
	for i := range err {
		if err[i] != nil {
			//Here I should add the line number to the returned error.
			return nil, 0, err[i]
		}
	}
	return coords, bfactor, nil
}

//PDBReadString reads the atomic entries for a PDB bufio.IO, returns a bunch of without coordinates,
// and the coordinates in a separate array of arrays. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not
// really well set up right now.
//WARNING: It has not been tested with an actual string.
func PDBStringRead(pdb string, read_additional bool) (*Molecule, error) {
	pdbstringreader := strings.NewReader(pdb)
	bufiopdb := bufio.NewReader(pdbstringreader)
	mol, err := pdbBufIORead(bufiopdb, read_additional)
	return mol, err
}

func PDBReaderRead(pdb io.Reader, read_additional bool) (*Molecule, error) {
	bufiopdb := bufio.NewReader(pdb)
	mol, err := pdbBufIORead(bufiopdb, read_additional)
	return mol, err
}

//PDBRead reads the atomic entries for a PDB file, returns a bunch of without coordinates,
// and the coordinates in a separate array of arrays. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not
// really well set up right now.
func PDBRead(pdbname string, read_additional bool) (*Molecule, error) {
	pdbfile, err := os.Open(pdbname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, err
	}
	defer pdbfile.Close()
	pdb := bufio.NewReader(pdbfile)
	mol, err := pdbBufIORead(pdb, read_additional)
	return mol, err
}

//pdbBufIORead reads the atomic entries for a PDB bufio.IO, returns a bunch of without coordinates,
// and the coordinates in a separate array of arrays. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not
// really well set up right now.
func pdbBufIORead(pdb *bufio.Reader, read_additional bool) (*Molecule, error) {
	molecule := make([]*Atom, 0)
	modelnumber := 0 //This is the number of frames read
	coords := make([][]float64, 1, 1)
	coords[0] = make([]float64, 0)
	bfactors := make([][]float64, 1, 1)
	bfactors[0] = make([]float64, 0)
	first_model := true //are we reading the first model? if not we only save coordinates
	contlines := 1      //count the lines read to better report errors
	for {
		line, err := pdb.ReadString('\n')
		if err != nil {
			//fmt.Println("PDB reading complete") /***change this to stderr************/
			break
			contlines++ //count all the lines even if empty.
		}
		if len(line) < 4 {
			continue
		}
		//here we start actually reading
		/*There might be a bug for not copying the string (symbol, name, etc) but just assigning the slice
		 * which is a reference. Check!*/
		var c = make([]float64, 3, 3)
		var bfactemp float64 //temporary bfactor
		var atomtmp *Atom
		//	var foo string // not really needed
		if strings.HasPrefix(line, "ATOM") || strings.HasPrefix(line, "HETATM") {
			if !first_model {
				c, bfactemp, err = read_onlycoords_pdb_line(line, contlines)
			} else {
				atomtmp = new(Atom)
				atomtmp, c, bfactemp, err = read_full_pdb_line(line, read_additional, contlines)
				if err != nil {
					return nil, err
				}
				//atom data other than coords is the same in all models so just read for the first.
				molecule = append(molecule, atomtmp)
			}
			//coords are appended for all the models
			//we add the coords to the latest frame of coordinaates
			coords[len(coords)-1] = append(coords[len(coords)-1], c[0], c[1], c[2])
			bfactors[len(bfactors)-1] = append(bfactors[len(bfactors)-1], bfactemp)
		} else if strings.HasPrefix(line, "MODEL") {
			modelnumber++        //,_=strconv.Atoi(strings.TrimSpace(line[6:]))
			if modelnumber > 1 { //will be one for the first model, 2 for the second.
				first_model = false
				coords = append(coords, make([]float64, 0)) //new bunch of coords for a new frame
				bfactors = append(bfactors, make([]float64, 0))
			}
		}
	}
	//This could be done faster if done in the same loop where the coords are read
	//Instead of having another loop just for them.
	top, err := NewTopology(molecule, 0, 0)
	if err != nil {
		return nil, err
	}
	frames := len(coords)
	mcoords := make([]*VecMatrix, frames, frames) //Final thing to return
	for i := 0; i < frames; i++ {
		mcoords[i], err = NewVecs(coords[i])
		if err != nil {
			return nil, err
		}
	}
	//if something happened during the process
	if err != nil {
		return nil, err
	}
	returned, err := NewMolecule(mcoords, top, bfactors)
	return returned, err
}

//End PDB_read family

//correctBfactors check that coords and bfactors have the same number of elements.
func correctBfactors(coords []*VecMatrix, bfactors [][]float64) bool {
	if len(coords) != len(bfactors) || bfactors==nil{
		return false
	}
	for key, coord := range coords {
		cr, _ := coord.Dims()
		br := len(bfactors[key])
		if cr != br {
			return false
		}
	}
	return true
}

//writePDBLine writes a line in PDB format from the data passed as a parameters. It takes the chain of the previous atom
//and returns the written line, the chain of the just-written atom, and error or nil.
func writePDBLine(atom *Atom, coord *VecMatrix, bfact float64, chainprev string) (string, string, error) {
	var ter string
	var out string
	if atom.Chain != chainprev {
		ter = fmt.Sprint(out, "TER\n")
		chainprev = atom.Chain
	}
	first := "ATOM"
	if atom.Het {
		first = "HETATM"
	}
	formatstring := "%-6s%5d  %-3s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"
	//4 chars for the atom name are used when hydrogens are included.
	//This has not been tested
	if len(atom.Name) == 4 {
		formatstring = strings.Replace(formatstring, "%-3s ", "%-4s", 1)
	} else if len(atom.Name) > 4 {
		return "", chainprev, fmt.Errorf("Cant print PDB line")
	}
	//"%-6s%5d  %-3s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"
	out = fmt.Sprintf(formatstring, first, atom.Id, atom.Name, atom.Molname, atom.Chain,
		atom.Molid, coord.At(0, 0), coord.At(0, 1), coord.At(0, 2), atom.Occupancy, bfact, atom.Symbol)

	out = strings.Join([]string{ter, out}, "")
	return out, chainprev, nil
}

//PDBWrite writes a PDB for the molecule mol and the coordinates Coords. It is just a wrapper for
//PDBStringWrite. Returns error or nil.
func PDBWrite(pdbname string, coords *VecMatrix, mol Atomer, Bfactors []float64) error {
	out, err := os.Create(pdbname)
	if err != nil {
		return err
	}
	defer out.Close()
	fmt.Fprint(out, "REMARK     WRITTEN WITH GOCHEM :-)\n")
	outstring, err := PDBStringWrite(coords, mol, Bfactors)
	if err != nil {
		return err
	}
	_, err = fmt.Fprintf(out, outstring) //this might require checking for the error
	if err != nil {
		return err
	}

	return nil
}

//PDBStringWrite writes a string in PDB format for a given reference, coordinate set and bfactor set, which must match each other
//returns the written string and error or nil.
func PDBStringWrite(coords *VecMatrix, mol Atomer, bfact []float64) (string, error) {
	if bfact == nil {
		bfact = make([]float64, mol.Len())
	}
	cr, _ := coords.Dims()
	br := len(bfact)
	if cr != mol.Len() || cr != br {
		return "", fmt.Errorf("Ref and Coords and/or Bfactors dont have the same number of atoms")
	}
	chainprev := mol.Atom(0).Chain //this is to know when the chain changes.
	var outline string
	var outstring string
	var err error
	for i := 0; i < mol.Len(); i++ {
		//	r,c:=coords.Dims()
		//	fmt.Println("IIIIIIIIIIIi", i,coords,r,c, "lllllll")
		writecoord := coords.VecView(i)
		outline, chainprev, err = writePDBLine(mol.Atom(i), writecoord, bfact[i], chainprev)
		if err != nil {
			return "", fmt.Errorf("Could not print PDB line: %d", i)
		}
		outstring = strings.Join([]string{outstring, outline}, "")
	}
	outstring = strings.Join([]string{outstring, "END\n"}, "")
	return outstring, nil
}

//MultiPDBWrite writes a multiPDB file for the molecule mol and the various coordinate sets in CandB.
//CandB is a list of lists of *matrix.DenseMatrix. If it has 2 elements or more, the second will be used as
//Bfactors. If it has one element, all b-factors will be zero.
//Returns an error if fails, or nil if succeeds.
func MultiPDBWrite(pdbname string, Coords []*VecMatrix, mol Atomer, Bfactors [][]float64) error {
	if !correctBfactors(Coords, Bfactors) {
		Bfactors = make([][]float64, len(Coords), len(Coords))
	}

	out, err := os.Create(pdbname)
	if err != nil {
		return err
	}
	defer out.Close()

	fmt.Fprint(out, "REMARK     WRITTEN WITH GOCHEM :-)\n")
	for j := range Coords {
		fmt.Fprintf(out, "MODEL %d\n", j+1) //The model number starts with one
		outstring, err := PDBStringWrite(Coords[j], mol, Bfactors[j])
		if err != nil {
			return err
		}
		outstring=strings.Replace(outstring, "END\n", "ENDMDL\n", 1)
		fmt.Fprintf(out, outstring) //here it could be needed to check for the error
	}

	fmt.Fprint(out, "END\n")
	return nil
}

/***End of PDB part***/

//Reads a string formated as an xyz or multixyz (as produced by Turbomole). Returns a Molecule and error or nil.
func XYZStringRead(xyz string) (*Molecule, error) {
	xyzstringreader := strings.NewReader(xyz)
	bufioxyz := bufio.NewReader(xyzstringreader)
	mol, err := XYZBufIORead(bufioxyz)
	return mol, err
}

//Reads an xyz or multixyz file (as produced by Turbomole). Returns a Molecule and error or nil.
func XYZRead(xyzname string) (*Molecule, error) {
	xyzfile, err := os.Open(xyzname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, err
	}
	defer xyzfile.Close()
	xyz := bufio.NewReader(xyzfile)
	mol, err := XYZBufIORead(xyz)
	if err != nil {
		errstr := err.Error()
		err = fmt.Errorf(strings.Join([]string{errstr, " in file ", xyzname}, ""))
	}
	return mol, err

}

//Reads an xyz or multixyz formatted bufio.Reader (as produced by Turbomole). Returns a Molecule and error or nil.
func XYZBufIORead(xyz *bufio.Reader) (*Molecule, error) {
	snaps := 1
	var err error
	var top *Topology
	var molecule []*Atom
	Coords := make([]*VecMatrix, 1, 1)

	for {
		//When we read the first snapshot we collect also the topology data, later
		//only coords are collected.
		if snaps == 1 {
			Coords[0], molecule, err = xyzReadSnap(xyz, true)
			if err != nil {
				return nil, err
			}
			top, err = NewTopology(molecule, 0, 0)
			if err != nil {
				return nil, err
			}
			snaps++
			continue
		}
		tmpcoords, _, err := xyzReadSnap(xyz, false)
		if err != nil {
			//An error here simply means that there are no more snapshots
			errm := err.Error()
			if strings.Contains(errm, "Empty") || strings.Contains(errm, "header") {
				err = nil
				break
			}
			return nil, err
		}
		Coords = append(Coords, tmpcoords)
	}
	bfactors := make([][]float64, len(Coords), len(Coords))
	for key, _ := range bfactors {
		bfactors[key] = make([]float64, top.Len())
	}
	returned, err := NewMolecule(Coords, top, bfactors)
	return returned, err
}

//xyzReadSnap reads an xyz file snapshot from a bufio.Reader, returns a slice of Atom objects, which will be nil if ReadTopol is false,
// a slice of matrix.DenseMatrix and an error or nil.
func xyzReadSnap(xyz *bufio.Reader, ReadTopol bool) (*VecMatrix, []*Atom, error) {
	line, err := xyz.ReadString('\n')
	if err != nil {
		return nil, nil, fmt.Errorf("Empty XYZ File: %s", err.Error())
	}
	natoms, err := strconv.Atoi(strings.TrimSpace(line))
	if err != nil {
		return nil, nil, fmt.Errorf("Wrong header for an XYZ file %s", err.Error())
	}
	var molecule []*Atom
	if ReadTopol {
		molecule = make([]*Atom, natoms, natoms)
	}
	coords := make([]float64, natoms*3, natoms*3)
	_, err = xyz.ReadString('\n') //We dont care about this line
	if err != nil {
		return nil, nil, fmt.Errorf("Ill formatted XYZ file: %s", err.Error())
	}
	errs := make([]error, 3, 3)
	for i := 0; i < natoms; i++ {
		line, errs[0] = xyz.ReadString('\n')
		if errs[0] != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			if strings.Contains(errs[0].Error(), "EOF") && i == natoms-1 { //This allows that an XYZ ends without a newline
				errs[0] = nil
			} else {
				break
			}
		}
		fields := strings.Fields(line)
		if len(fields) < 4 {
			errs[0] = fmt.Errorf("Line number %d ill formed", i)
			break
		}
		if ReadTopol {
			molecule[i] = new(Atom)
			molecule[i].Symbol = strings.Title(fields[0])
			molecule[i].Mass = symbolMass[molecule[i].Symbol]
			molecule[i].Molname = "UNK"
			molecule[i].Name = molecule[i].Symbol
		}
		coords[i*3], errs[0] = strconv.ParseFloat(fields[1], 64)
		coords[i*3+1], errs[1] = strconv.ParseFloat(fields[2], 64)
		coords[i*3+2], errs[2] = strconv.ParseFloat(fields[3], 64)
	}
	//This could be done faster if done in the same loop where the coords are read
	//Instead of having another loop just for them.
	for k, i := range errs {
		if i != nil {
			fmt.Println("line", line, k)
			return nil, nil, i
		}
	}
	mcoords, err := NewVecs(coords)
	return mcoords, molecule, err
}

//XYZWrite writes the mol Ref and the Coord coordinates in an XYZ file with name xyzname which will
//be created fot that. If the file exist it will be overwritten.
func XYZWrite(xyzname string, Coords *VecMatrix, mol Atomer) error {
	out, err := os.Create(xyzname)
	if err != nil {
		return err
	}
	defer out.Close()
	xyz, err := XYZStringWrite(Coords, mol)
	if err != nil {
		return err
	}
	fmt.Fprintf(out, xyz)
	return nil
}

//XYZStringWrite writes the mol Ref and the Coord coordinates in an XYZ-formatted string.
func XYZStringWrite(Coords *VecMatrix, mol Atomer) (string, error) {
	var out string
	if mol.Len() != Coords.NVecs() {
		return "", fmt.Errorf("Ref and Coords dont have the same number of atoms")
	}
	c := make([]float64, 3, 3)
	out = fmt.Sprintf("%-4d\n\n", mol.Len())
	//towrite := Coords.Arrays() //An array of array with the data in the matrix
	for i := 0; i < mol.Len(); i++ {
		//c := towrite[i] //coordinates for the corresponding atoms
		c = Coords.Row(c, i)
		temp := fmt.Sprintf("%-2s  %12.6f%12.6f%12.6f \n", mol.Atom(i).Symbol, c[0], c[1], c[2])
		out = strings.Join([]string{out, temp}, "")
	}
	return out, nil
}
