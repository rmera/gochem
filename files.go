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

import "os"
import "bufio"
import "strings"
import "strconv"
import "fmt"
import "github.com/skelterjohn/go.matrix"

//Pdb_read family

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
	atom.Chain = line[21] //currently this is read to a byte
	atom.Molid, err[1] = strconv.Atoi(strings.TrimSpace(line[22:26]))
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
		atom.Symbol = strings.Title(strings.ToLower(atom.Symbol)) //Not too efficient I guess
		atom.Charge = float64(line[78])                           //strconv.ParseFloat(strings.TrimSpace(line[78:78]),64)
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

//PdbRead reads the atomic entries for a PDB file, returns a bunch of without coordinates,
// and the coordinates in a separate array of arrays. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not 
// really well set up right now.
func PdbRead(pdbname string, read_additional bool) (*Molecule, error) {
	molecule := make([]*Atom, 0)
	modelnumber := 0 //This is the number of frames read
	coords := make([][]float64, 1, 1)
	coords[0] = make([]float64, 0)
	bfactors := make([][]float64, 1, 1)
	bfactors[0] = make([]float64, 0)
	first_model := true //are we reading the first model? if not we only save coordinates
	pdbfile, err := os.Open(pdbname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, err
	}
	defer pdbfile.Close()
	pdb := bufio.NewReader(pdbfile)
	contlines := 1 //count the lines read to better report errors
	for {
		line, err := pdb.ReadString('\n')
		if err != nil {
			fmt.Println("PDB reading complete") /***change this to stderr************/
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
	top, err := MakeTopology(molecule, 0, 0)
	if err != nil {
		return nil, err
	}
	frames := len(coords)
	mcoords := make([]*matrix.DenseMatrix, frames, frames) //Final thing to return
	mbfactors := make([]*matrix.DenseMatrix, frames, frames)
	for i := 0; i < frames; i++ {
		mcoords[i] = matrix.MakeDenseMatrix(coords[i], len(coords[i])/3, 3)
		mbfactors[i] = matrix.MakeDenseMatrix(bfactors[i], len(bfactors[i]), 1)
	}
	//if something happened during the process
	if err != nil {
		return nil, err
	}
	returned, err := MakeMolecule(top, mcoords, mbfactors)
	return returned, err
}

//End Pdb_read family

//PdbWrite writes a PDB for the molecule mol and the coordinates Coords. It is just a wrapper for
//MultiPdbWrite. Returns error or nil.
func PdbWrite(pdbname string, mol Ref, CandB ...*matrix.DenseMatrix) error {
	Coords := []*matrix.DenseMatrix{CandB[0]}
	var Bfactors []*matrix.DenseMatrix
	if len(CandB) > 1 {
		Bfactors = []*matrix.DenseMatrix{CandB[1]} //any other element is just ignored	
	} else {
		Bfactors = []*matrix.DenseMatrix{matrix.Zeros(mol.Len(), 1)}
	}
	err := MultiPdbWrite(pdbname, mol, Coords, Bfactors)
	return err
}

//MultiPdbWrite writes a multiPDB file for the molecule mol and the various coordinate sets in CandB.
//CandB is a list of lists of *matrix.DenseMatrix. If it has 2 elements or more, the second will be used as
//Bfactors. If it has one element, all b-factors will be zero.
//Returns an error if fails, or nil if succeeds.
func MultiPdbWrite(pdbname string, mol Ref, CandB ...[]*matrix.DenseMatrix) error {
	/*if len(CandB)<1{  //We should not need this
		return fmt.Errorf("MultiPdbWrite: Coordinates not supplied")
	}*/
	Coords := CandB[0]
	var Bfactors []*matrix.DenseMatrix
	if len(CandB) > 1 {
		Bfactors = CandB[1] //any other element is just ignored	
	} else {
		for _, _ = range Coords {
			Bfactors = append(Bfactors, matrix.Zeros(mol.Len(), 1))
		}
	}

	for key, val := range Coords {
		if val.Rows() != mol.Len() || val.Rows() != Bfactors[key].Rows() {
			return fmt.Errorf("MultiPdbWrite: Ref and Coords and/or Bfactors dont have the same number of atoms")
		}
	}
	out, err := os.Create(pdbname)
	if err != nil {
		return err
	}

	defer out.Close()
	fmt.Fprint(out, "REMARK     WRITTEN WITH GOCHEM :-)\n")
	for j := range Coords {
		towrite := Coords[j].Arrays()       //An array of array with the data in the matrix
		chainprev := mol.Atom(0).Chain      //this is to know when the chain changes.
		fmt.Fprintf(out, "MODEL %d\n", j+1) //The model number starts with one
		//		fmt.Println("chain", mol.Atoms[0])		
		for i := 0; i < mol.Len(); i++ {
			ThisAtom := mol.Atom(i)
			if ThisAtom.Chain != chainprev {
				fmt.Fprintln(out, "TER")
				chainprev = ThisAtom.Chain
			}
			first := "ATOM"
			if ThisAtom.Het {
				first = "HETATM"
			}
			//Corrupted will check for broken bfactors and complete with zeroes
			//So no need to check here.
			c := towrite[i] //coordinates for the corresponding atoms
			if len(ThisAtom.Name) < 4 {
				//		fmt.Println("chain", ThisAtom)
				_, err = fmt.Fprintf(out, "%-6s%5d  %-3s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n", first, ThisAtom.Id, ThisAtom.Name, ThisAtom.Molname, ThisAtom.Chain,
					ThisAtom.Molid, c[0], c[1], c[2], ThisAtom.Occupancy, Bfactors[j].Get(i, 0), ThisAtom.Symbol)
				//4 chars for the atom name are used when hydrogens are included.	
				//This has not been tested
			} else if len(ThisAtom.Name) == 4 {
				_, err = fmt.Fprintf(out, "%-6s%5d %4s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n", first, ThisAtom.Id, ThisAtom.Name, ThisAtom.Molname, ThisAtom.Chain,
					ThisAtom.Molid, c[0], c[1], c[2], ThisAtom.Occupancy, Bfactors[j].Get(i, 0), ThisAtom.Symbol)
			} else {
				err = fmt.Errorf("Cant print PDB line")
			}
			if err != nil {
				return err
			}
		}
		fmt.Fprint(out, "ENDMDL\n")
	}
	fmt.Fprint(out, "END\n")
	return nil
}

//Reads an xyz or multixyz file (as produced by Turbomole). Returns a Molecule and error or nil.
func XyzRead(xyzname string) (*Molecule, error) {
	xyzfile, err := os.Open(xyzname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, err
	}
	defer xyzfile.Close()
	Coords := make([]*matrix.DenseMatrix, 1, 1)
	xyz := bufio.NewReader(xyzfile)
	snaps := 1
	var top *Topology
	var molecule []*Atom
	for {
		//When we read the first snapshot we collect also the topology data, later
		//only coords are collected.
		if snaps == 1 {
			Coords[0], molecule, err = xyzReadSnap(xyz, xyzname, true)
			if err != nil {
				return nil, err
			}
			top, err = MakeTopology(molecule, 0, 0)
			if err != nil {
				return nil, err
			}
			snaps++
			continue
		}
		tmpcoords, _, err := xyzReadSnap(xyz, xyzname, false)
		if err != nil {
			//An error here simply means that there are no more snapshots
			err = nil
			break
		}
		Coords = append(Coords, tmpcoords)
	}
	bfactors := make([]*matrix.DenseMatrix, len(Coords), len(Coords))
	for key, _ := range bfactors {
		bfactors[key] = matrix.Zeros(top.Len(), 1)
	}
	returned, err := MakeMolecule(top, Coords, bfactors)
	return returned, err
}

//XyzRead reads an xyz file, returns a slice of Atom objects, which will be nil if ReadTopol is false,
// a slice of matrix.DenseMatrix and an error or nil.
func xyzReadSnap(xyz *bufio.Reader, xyzname string, ReadTopol bool) (*matrix.DenseMatrix, []*Atom, error) {
	line, err := xyz.ReadString('\n')
	if err != nil {
		return nil, nil, fmt.Errorf("Ill formatted XYZ file")
	}
	natoms, err := strconv.Atoi(strings.TrimSpace(line))
	if err != nil {
		return nil, nil, fmt.Errorf("Ill formatted XYZ file")
	}
	var molecule []*Atom
	if ReadTopol {
		molecule = make([]*Atom, natoms, natoms)
	}
	coords := make([]float64, natoms*3, natoms*3)
	_, err = xyz.ReadString('\n') //We dont care about this line
	if err != nil {
		return nil, nil, fmt.Errorf("Ill formatted XYZ file")
	}
	errs := make([]error, 3, 3)
	for i := 0; i < natoms; i++ {
		line, errs[0] = xyz.ReadString('\n')
		if errs[0] != nil { //inefficient, (errs[1] can be checked once before), but clearer.
			break
		}
		fields := strings.Fields(line)
		if len(fields) < 4 {
			errs[0] = fmt.Errorf("Line number %d in file %s ill formed", i, xyzname)
			break
		}
		if ReadTopol {
			molecule[i] = new(Atom)
			molecule[i].Symbol = strings.ToTitle(strings.ToLower(fields[0]))
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
	for _, i := range errs {
		if i != nil {
			return nil, nil, i
		}
	}
	mcoords := matrix.MakeDenseMatrix(coords, natoms, 3)
	return mcoords, molecule, err
}

//XyzWrite writes the mol Ref and the Coord coordinates in an XYZ file with name xyzname which will
//be created fot that. If the file exist it will be overwritten.
func XyzWrite(xyzname string, mol Ref, Coords *matrix.DenseMatrix) error {
	if mol.Len() != Coords.Rows() {
		return fmt.Errorf("Ref and Coords dont have the same number of atoms")
	}
	out, err := os.Create(xyzname)
	if err != nil {
		return err
	}
	defer out.Close()
	fmt.Fprintf(out, "%-4d\n", mol.Len())
	fmt.Fprintf(out, "\n")
	towrite := Coords.Arrays() //An array of array with the data in the matrix	
	for i := 0; i < mol.Len(); i++ {
		c := towrite[i] //coordinates for the corresponding atoms
		_, err = fmt.Fprintf(out, "%-2s  %12.6f%12.6f%12.6f \n", mol.Atom(i).Symbol, c[0], c[1], c[2])
		if err != nil {
			return err
		}
	}
	return nil
}
