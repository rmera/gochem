/*
 * files.go, part of gochem.
 *
 *
 * Copyright 2012 Raul Mera rauldotmeraatusachdotcl
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
 * goChem is developed at Universidad de Santiago de Chile (USACH)
 *
 *
 */

package chem

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	v3 "github.com/rmera/gochem/v3"
)

//PDB_read family

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
	} else if name[0] == 'B' {
		if name == "BE" {
			symbol = "Be"
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
		return symbol, CError{"Couldn't guess symbol from PDB name", []string{"symbolFromName"}}
	}
	return symbol, nil
}

// read_full_pdb_line parses a valid ATOM or HETATM line of a PDB file, returns an Atom
// object with the info except for the coordinates and b-factors, which  are returned
// separately as an array of 3 float64 and a float64, respectively
func read_full_pdb_line(line string, read_additional bool, contlines int) (*Atom, []float64, float64, error) {
	err := make([]error, 7, 7) //accumulate errors to check at the end of the readed line.
	coords := make([]float64, 3, 3)
	atom := new(Atom)
	atom.Het = strings.HasPrefix(line, "HETATM") //this is called twice in the worst case. should fix
	atom.ID, err[0] = strconv.Atoi(strings.TrimSpace(line[6:12]))
	atom.Name = strings.TrimSpace(line[12:16])
	atom.Char16 = line[16]
	//PDB says that pos. 17 is for other thing but I see that is
	//used for residue name in many cases*/
	atom.MolName = line[17:20]
	atom.MolName1 = three2OneLetter[atom.MolName]
	atom.Chain = string(line[21])
	atom.MolID, err[1] = strconv.Atoi(strings.TrimSpace(line[22:30]))
	//Here we shouldn't need TrimSpace, but I keep it just in case someone
	// doesn's use all the fields when writting a PDB*/
	coords[0], err[2] = strconv.ParseFloat(strings.TrimSpace(line[30:38]), 64)
	coords[1], err[3] = strconv.ParseFloat(strings.TrimSpace(line[38:46]), 64)
	coords[2], err[4] = strconv.ParseFloat(strings.TrimSpace(line[46:54]), 64)
	var bfactor float64
	//Every correct PDB should include occupancy and b-factor, but _of course_ writing
	//correct PDBs is too hard for some programs (and by "some programs" I mean OPLS LigParGen. Get it together, guys).
	//so I add this conditional to allow goChem to still read these wrong PDB files.
	if len(line) >= 60 {
		atom.Occupancy, err[5] = strconv.ParseFloat(strings.TrimSpace(line[54:60]), 64)
		bfactor, err[6] = strconv.ParseFloat(strings.TrimSpace(line[60:66]), 64)
	}
	//we try to read the additional only if indicated and if it is there
	// In this part we don't catch errors. If something is missing we
	// just ommit it
	if read_additional && len(line) >= 79 {
		atom.Symbol = strings.TrimSpace(line[76:78])
		atom.Symbol = strings.Title(strings.ToLower(atom.Symbol))
		if len(line) >= 80 {
			var errcharge error
			atom.Charge, errcharge = strconv.ParseFloat(strings.TrimSpace(line[78:78]), 64)
			if errcharge == nil {

				if strings.Contains(line[79:79], "-") {
					atom.Charge = -1.0 * atom.Charge
				}
			} else {
				//we dont' report an error here, just set the charge to 0 (the default)
				atom.Charge = 0.0
			}
		}
	}

	//This part tries to guess the symbol from the atom name, if it has not been read
	//No error checking here, just fills symbol with the empty string the function returns
	var symbolerr error
	if len(atom.Symbol) == 0 {
		atom.Symbol, symbolerr = symbolFromName(atom.Name)
	}

	for i := range err {
		if err[i] != nil {
			//Here I should add the line number to the returned error.
			return nil, nil, 0, CError{err[i].Error(), []string{"strconv.Atoi/ParseFloat", "read_full_pdb_line"}}
		}
	}
	if atom.Symbol != "" {
		atom.Mass = symbolMass[atom.Symbol] //Not error checking
	}
	//if we couldn't read the symbol, we'll still return the atom and coords
	//but with an error
	return atom, coords, bfactor, symbolerr
}

//read_onlycoords_pdb_line parses an ATOM/HETATM PDB line returning only the coordinates and b-factors
func read_onlycoords_pdb_line(line string, contlines int) ([]float64, float64, error) {
	coords := make([]float64, 3, 3)
	err := make([]error, 4, 4)
	var bfactor float64
	coords[0], err[0] = strconv.ParseFloat(strings.TrimSpace(line[30:38]), 64)
	coords[1], err[1] = strconv.ParseFloat(strings.TrimSpace(line[38:46]), 64)
	coords[2], err[2] = strconv.ParseFloat(strings.TrimSpace(line[46:54]), 64)
	bfactor, err[3] = strconv.ParseFloat(strings.TrimSpace(line[60:66]), 64)
	for i := range err {
		if err[i] != nil {
			//Here I should add the line number to the returned error.
			return nil, 0, CError{err[i].Error(), []string{"strconv.ParseFloat", "read_onlycoords_pdb_line"}}

		}
	}
	return coords, bfactor, nil
}

//PDBRRead reads a pdb file from an io.Reader. Returns a Molecule. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not
// really well set up right now.
//read_additional is now "deprecated", it will be set to true regardless. I have made it into
func PDBRead(pdb io.Reader, read_additional ...bool) (*Molecule, error) {
	bufiopdb := bufio.NewReader(pdb)
	mol, err := pdbBufIORead(bufiopdb, read_additional...)
	return mol, errDecorate(err, "PDBReaderREad")
}

//PDBFileRead reads a pdb file from an io.Reader. Returns a Molecule. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not
// really well set up right now. read_additional is now deprecated. The reader will just read
func PDBFileRead(pdbname string, read_additional ...bool) (*Molecule, error) {
	pdbfile, err := os.Open(pdbname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, err
	}
	defer pdbfile.Close()
	pdb := bufio.NewReader(pdbfile)
	mol, err := pdbBufIORead(pdb, read_additional...)
	return mol, err
}

//pdbBufIORead reads the atomic entries for a PDB bufio.IO, reads a pdb file from an io.Reader.
//Returns a Molecule. If there is one frame in the PDB the coordinates array will be of lenght 1.
//It also returns an error which is not really well set up right now.
//The read_additional_opt allows not reading the last fields of a PDB, if you know they are wrong.
//if true (the default), the fields are read if they are available. Otherwise we attempt to figure
//out the symbol from the atom name, which doesn't always work.
func pdbBufIORead(pdb *bufio.Reader, read_additional_opt ...bool) (*Molecule, error) {
	read_additional := true
	if len(read_additional_opt) > 0 {
		read_additional = read_additional_opt[0]
	}
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
			//	contlines++ //count all the lines even if empty. This is unreachable but I'm not sure at this point if it's better this way! goChem does read PDBs correctly as far as I can see.
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
				if err != nil {
					return nil, errDecorate(err, "pdbBufIORead")
				}
			} else {
				atomtmp = new(Atom)
				atomtmp, c, bfactemp, err = read_full_pdb_line(line, read_additional, contlines)
				if err != nil {
					return nil, errDecorate(err, "pdbBufIORead")
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
	top := NewTopology(0, 1, molecule)
	var err error
	frames := len(coords)
	mcoords := make([]*v3.Matrix, frames, frames) //Final thing to return
	for i := 0; i < frames; i++ {
		mcoords[i], err = v3.NewMatrix(coords[i])
		if err != nil {
			return nil, errDecorate(err, "pdbBufIORead")
		}
	}
	//if something happened during the process
	if err != nil {
		return nil, errDecorate(err, "pdbBufIORead")
	}
	returned, err := NewMolecule(mcoords, top, bfactors)
	return returned, errDecorate(err, "pdbBufIORead")
}

//End PDB_read family

//correctBfactors check that coords and bfactors have the same number of elements.
func correctBfactors(coords []*v3.Matrix, bfactors [][]float64) bool {
	if len(coords) != len(bfactors) || bfactors == nil {
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
func writePDBLine(atom *Atom, coord *v3.Matrix, bfact float64, chainprev string) (string, string, error) {
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
		return "", chainprev, CError{"Cant print PDB line", []string{"writePDBLine"}}
	}
	//"%-6s%5d  %-3s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"
	out = fmt.Sprintf(formatstring, first, atom.ID, atom.Name, atom.MolName, atom.Chain,
		atom.MolID, coord.At(0, 0), coord.At(0, 1), coord.At(0, 2), atom.Occupancy, bfact, atom.Symbol)
	out = strings.Join([]string{ter, out}, "")
	return out, chainprev, nil
}

//PDBFileWrite writes a PDB for the molecule mol and the coordinates Coords to a file name pdbname.
func PDBFileWrite(pdbname string, coords *v3.Matrix, mol Atomer, Bfactors []float64) error {
	out, err := os.Create(pdbname)
	if err != nil {
		return CError{err.Error(), []string{"os.Create", "PDBFileWrite"}}
	}
	defer out.Close()
	fmt.Fprintf(out, "REMARK WRITTEN WITH GOCHEM :-) \n")
	err = PDBWrite(out, coords, mol, Bfactors)
	if err != nil {
		return errDecorate(err, "PDBFileWrite")
	}
	return nil
}

//PDBWrite writes a PDB formatted sequence of bytes to an io.Writer for a given reference, coordinate set and bfactor set, which must match each other. Returns error or nil.
func PDBWrite(out io.Writer, coords *v3.Matrix, mol Atomer, bfact []float64) error {
	err := pdbWrite(out, coords, mol, bfact)
	if err != nil {
		errDecorate(err, "PDBWrite")
	}
	_, err = out.Write([]byte{'\n'}) //This function is just a wrapper to add the newline to what pdbWrite does.
	if err != nil {
		return CError{"Failed to write in io.Writer", []string{"io.Write.Write", "PDBWrite"}}

	}
	return nil
}

func pdbWrite(out io.Writer, coords *v3.Matrix, mol Atomer, bfact []float64) error {
	if bfact == nil {
		bfact = make([]float64, mol.Len())
	}
	cr, _ := coords.Dims()
	br := len(bfact)
	if cr != mol.Len() || cr != br {
		return CError{"Ref and Coords and/or Bfactors dont have the same number of atoms", []string{"pdbWrite"}}
	}
	chainprev := mol.Atom(0).Chain //this is to know when the chain changes.
	var outline string
	var err error
	iowriteError := func(err error) error {
		return CError{"Failed to write in io.Writer" + err.Error(), []string{"io.Write.Write", "pdbWrite"}}
	}
	for i := 0; i < mol.Len(); i++ {
		//	r,c:=coords.Dims()
		//	fmt.Println("IIIIIIIIIIIi", i,coords,r,c, "lllllll")
		writecoord := coords.VecView(i)
		outline, chainprev, err = writePDBLine(mol.Atom(i), writecoord, bfact[i], chainprev)
		if err != nil {
			return errDecorate(err, "pdbWrite "+fmt.Sprintf("Could not print PDB line: %d", i))
		}
		_, err := out.Write([]byte(outline))
		if err != nil {
			return iowriteError(err)
		}
	}
	_, err = out.Write([]byte("TER\n")) // New Addition, should help to recognize the end of the chain.
	_, err = out.Write([]byte("END"))   //no newline, this is in case the write is part of a PDB and one needs to write "ENDMDEL".
	if err != nil {
		return iowriteError(err)
	}
	return nil
}

//PDBStringWrite writes a string in PDB format for a given reference, coordinate set and bfactor set, which must match each other
//returns the written string and error or nil.
func PDBStringWrite(coords *v3.Matrix, mol Atomer, bfact []float64) (string, error) {
	if bfact == nil {
		bfact = make([]float64, mol.Len())
	}
	cr, _ := coords.Dims()
	br := len(bfact)
	if cr != mol.Len() || cr != br {
		return "", CError{"Ref and Coords and/or Bfactors dont have the same number of atoms", []string{"PDBStringWrite"}}
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
			return "", errDecorate(err, "PDBStringWrite "+fmt.Sprintf("Could not print PDB line: %d", i))
		}
		outstring = strings.Join([]string{outstring, outline}, "")
	}
	outstring = strings.Join([]string{outstring, "END\n"}, "")
	return outstring, nil
}

//MultiPDBWrite writes a multiPDB  for the molecule mol and the various coordinate sets in CandB, to the io.Writer given.
//CandB is a list of lists of *matrix.DenseMatrix. If it has 2 elements or more, the second will be used as
//Bfactors. If it has one element, all b-factors will be zero.
//Returns an error if fails, or nil if succeeds.
func MultiPDBWrite(out io.Writer, Coords []*v3.Matrix, mol Atomer, Bfactors [][]float64) error {
	if !correctBfactors(Coords, Bfactors) {
		Bfactors = make([][]float64, len(Coords), len(Coords))
	}
	iowriterError := func(err error) error {
		return CError{"Failed to write in io.Writer" + err.Error(), []string{"io.Writer.Write", "MultiPDBWrite"}}
	}

	_, err := out.Write([]byte("REMARK WRITTEN WITH GOCHEM :-)")) //The model number starts with one
	if err != nil {
		return iowriterError(err)
	}
	//OK now the real business.
	for j := range Coords {
		_, err := out.Write([]byte(fmt.Sprintf("MODEL %d\n", j+1))) //The model number starts with one
		if err != nil {
			return iowriterError(err)
		}
		err = pdbWrite(out, Coords[j], mol, Bfactors[j])
		if err != nil {
			return errDecorate(err, "MultiPDBWrite")
		}
		_, err = out.Write([]byte("MDL\n"))
		if err != nil {
			return iowriterError(err)
		}

	}

	_, err = out.Write([]byte("END\n"))
	if err != nil {
		return iowriterError(err)
	}

	return nil
}

/***End of PDB part***/

//XYZFileRead Reads an xyz or multixyz file (as produced by Turbomole). Returns a Molecule and error or nil.
func XYZFileRead(xyzname string) (*Molecule, error) {
	xyzfile, err := os.Open(xyzname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, CError{err.Error(), []string{"os.Open", "XYZFileRead"}}
	}
	defer xyzfile.Close()
	mol, err := XYZRead(xyzfile)
	if err != nil {
		err = errDecorate(err, "XYZFileRead "+fmt.Sprintf(strings.Join([]string{"error in file ", xyzname}, "")))
	}
	return mol, err

}

//XYZRead Reads an xyz or multixyz formatted bufio.Reader (as produced by Turbomole). Returns a Molecule and error or nil.
func XYZRead(xyzp io.Reader) (*Molecule, error) {
	snaps := 1
	xyz := bufio.NewReader(xyzp)
	var err error
	var top *Topology
	var molecule []*Atom
	Coords := make([]*v3.Matrix, 1, 1)
	Data := make([]string, 1, 1)

	for {
		//When we read the first snapshot we collect also the topology data, later
		//only coords are collected.
		if snaps == 1 {
			Coords[0], molecule, Data[0], err = xyzReadSnap(xyz, nil, true)
			if err != nil {
				return nil, errDecorate(err, "XYZRead")
			}
			top = NewTopology(0, 1, molecule)
			if err != nil {
				return nil, errDecorate(err, "XYZRead")
			}
			snaps++
			continue
		}
		tmpcoords, _, data, err := xyzReadSnap(xyz, nil, false)
		if err != nil {
			//An error here simply means that there are no more snapshots
			errm := err.Error()
			if strings.Contains(errm, "Empty") || strings.Contains(errm, "header") {
				err = nil
				break
			}
			return nil, errDecorate(err, "XYZRead")
		}
		Coords = append(Coords, tmpcoords)
		Data = append(Data, data)
	}
	bfactors := make([][]float64, len(Coords), len(Coords))
	for key, _ := range bfactors {
		bfactors[key] = make([]float64, top.Len())
	}
	returned, err := NewMolecule(Coords, top, bfactors)
	returned.XYZFileData = Data
	return returned, errDecorate(err, "XYZRead")
}

//XYZTraj is a trajectory-like representation of an XYZ File.
type XYZTraj struct {
	natoms     int
	xyz        *bufio.Reader //The DCD file
	frames     int
	xyzfile    *os.File
	readable   bool
	firstframe *v3.Matrix
}

//Readable returns true if the trajectory is fit to be read, false otherwise.
func (X *XYZTraj) Readable() bool {
	return X.readable
}

func (X *XYZTraj) Len() int {
	return X.natoms
}

//xyztrajerror returns a LastFrameError if the given error message contains certain keywords.
//otherwise, returns the original error.
func (X *XYZTraj) xyztrajerror(err error) error {
	errm := err.Error()
	X.xyzfile.Close()
	X.readable = false
	if strings.Contains(errm, "Empty") || strings.Contains(errm, "header") {
		return newlastFrameError("", X.frames)
	} else {
		return err
	}

}

//Next reads the next snapshot of the trajectory into coords, or discards it, if coords
//is nil. It can take a box slice of floats, but won't do anything with it
//(only for compatibility with the Traj interface.
func (X *XYZTraj) Next(coords *v3.Matrix, box ...[]float64) error {
	if coords == nil {
		_, _, _, err := xyzReadSnap(X.xyz, coords, false)
		if err != nil {
			//An error here probably means that there are no more snapshots
			return X.xyztrajerror(err)
		}
		X.frames++
		return nil
	}
	if X.frames == 0 {
		coords.Copy(X.firstframe) //slow, but I don't want to mess with the pointer I got.
		X.frames++
		X.firstframe = nil
		return nil
	}
	_, _, _, err := xyzReadSnap(X.xyz, coords, false)
	if err != nil {
		//An error here probably means that there are no more snapshots
		return X.xyztrajerror(err)
	}

	X.frames++
	return err
}

//Reads a multi-xyz file. Returns the first snapshot as a molecule, and the other ones as a XYZTraj
func XYZFileAsTraj(xyzname string) (*Molecule, *XYZTraj, error) {
	xyzfile, err := os.Open(xyzname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, nil, CError{err.Error(), []string{"os.Open", "XYZFileRead"}}
	}
	xyz := bufio.NewReader(xyzfile)
	//the molecule first
	coords, atoms, _, err := xyzReadSnap(xyz, nil, true)
	top := NewTopology(0, 1, atoms)
	bfactors := make([][]float64, 1, 1)
	bfactors[0] = make([]float64, top.Len())
	returned, err := NewMolecule([]*v3.Matrix{coords}, top, bfactors)
	//now the traj
	traj := new(XYZTraj)
	traj.xyzfile = xyzfile
	traj.xyz = xyz
	traj.natoms = returned.Len()
	traj.readable = true
	traj.firstframe = coords
	return returned, traj, nil
}

//xyzReadSnap reads an xyz file snapshot from a bufio.Reader, returns a slice of Atom
//objects, which will be nil if ReadTopol is false,
// a slice of matrix.DenseMatrix and an error or nil.
func xyzReadSnap(xyz *bufio.Reader, toplace *v3.Matrix, ReadTopol bool) (*v3.Matrix, []*Atom, string, error) {
	line, err := xyz.ReadString('\n')
	if err != nil {
		return nil, nil, "", CError{fmt.Sprintf("Empty XYZ File: %s", err.Error()), []string{"bufio.Reader.ReadString", "xyzReadSnap"}}
	}
	natoms, err := strconv.Atoi(strings.TrimSpace(line))
	if err != nil {
		return nil, nil, "", CError{fmt.Sprintf("Wrong header for an XYZ file %s", err.Error()), []string{"strconv.Atoi", "xyzReadSnap"}}
	}
	var molecule []*Atom
	if ReadTopol {
		molecule = make([]*Atom, natoms, natoms)
	}
	var coords []float64
	if toplace == nil {
		coords = make([]float64, natoms*3, natoms*3)
	} else {
		coords = toplace.RawSlice()
	}
	data, err := xyz.ReadString('\n') //The text in "data" could be anything, including just "\n"
	if err != nil {
		return nil, nil, "", CError{fmt.Sprintf("Ill formatted XYZ file: %s", err.Error()), []string{"bufio.Reader.ReadString", "xyzReadSnap"}}

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
			errs[0] = fmt.Errorf("Line number %d ill formed: %s", i, line)
			break
		}
		if ReadTopol {
			molecule[i] = new(Atom)
			molecule[i].Symbol = strings.Title(fields[0])
			molecule[i].Mass = symbolMass[molecule[i].Symbol]
			molecule[i].MolName = "UNK"
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
			//	fmt.Println("line", line, k)
			return nil, nil, "", CError{i.Error(), []string{"strconv.ParseFloat", "xyzReadSnap"}}
		}
	}
	//this should be fine even if I had a toplace matrix. Both toplace and mcoord should just point to the same data.
	mcoords, err := v3.NewMatrix(coords)
	return mcoords, molecule, data, errDecorate(err, "xyzReadSnap")
}

//XYZWrite writes the mol Ref and the Coord coordinates in an XYZ file with name xyzname which will
//be created fot that. If the file exist it will be overwritten.
func XYZFileWrite(xyzname string, Coords *v3.Matrix, mol Atomer) error {
	out, err := os.Create(xyzname)
	if err != nil {
		return CError{err.Error(), []string{"os.Create", "XYZFileWrite"}}
	}
	defer out.Close()
	err = XYZWrite(out, Coords, mol)
	if err != nil {
		return errDecorate(err, "XYZFileWrite")
	}
	return nil
}

//XYZStringWrite writes the mol Ref and the Coord coordinates in an XYZ-formatted string.
func XYZStringWrite(Coords *v3.Matrix, mol Atomer) (string, error) {
	var out string
	if mol.Len() != Coords.NVecs() {
		return "", CError{"Ref and Coords dont have the same number of atoms", []string{"XYZStringWrite"}}
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

//XYZWrite writes the mol Ref and the Coords coordinates to a io.Writer, in the XYZ format.
func XYZWrite(out io.Writer, Coords *v3.Matrix, mol Atomer) error {
	iowriterError := func(err error) error {
		return CError{"Failed to write in io.Writer" + err.Error(), []string{"io.Writer.Write", "XYZWrite"}}
	}
	if mol.Len() != Coords.NVecs() {
		return CError{"Ref and Coords dont have the same number of atoms", []string{"XYZWrite"}}
	}
	c := make([]float64, 3, 3)
	_, err := out.Write([]byte(fmt.Sprintf("%-4d\n\n", mol.Len())))
	if err != nil {
		return iowriterError(err)
	}
	//towrite := Coords.Arrays() //An array of array with the data in the matrix
	for i := 0; i < mol.Len(); i++ {
		//c := towrite[i] //coordinates for the corresponding atoms
		c = Coords.Row(c, i)
		temp := fmt.Sprintf("%-2s  %12.6f%12.6f%12.6f \n", mol.Atom(i).Symbol, c[0], c[1], c[2])
		_, err := out.Write([]byte(temp))
		if err != nil {
			return iowriterError(err)
		}
	}
	return nil
}

//GroFileRead reads a file in the Gromacs gro format, returning a molecule.
func GroFileRead(groname string) (*Molecule, error) {
	grofile, err := os.Open(groname)
	if err != nil {
		//fmt.Println("Unable to open file!!")
		return nil, CError{err.Error(), []string{"os.Open", "GroFileRead"}}
	}
	defer grofile.Close()
	snaps := 1
	gro := bufio.NewReader(grofile)
	var top *Topology
	var molecule []*Atom
	Coords := make([]*v3.Matrix, 1, 1)

	for {
		//When we read the first snapshot we collect also the topology data, later
		//only coords are collected.
		if snaps == 1 {
			Coords[0], molecule, err = groReadSnap(gro, true)
			if err != nil {
				return nil, errDecorate(err, "GroFileRead")
			}
			top = NewTopology(0, 1, molecule)
			if err != nil {
				return nil, errDecorate(err, "GroFileRead")
			}
			snaps++
			continue
		}
		//fmt.Println("how manytimes?") /////////////////////
		tmpcoords, _, err := groReadSnap(gro, false)
		if err != nil {
			break //We just ignore errors after the first snapshot, and simply read as many snapshots as we can.
			/*
				//An error here may just mean that there are no more snapshots
				errm := err.Error()
				if strings.Contains(errm, "Empty") || strings.Contains(errm, "EOF") {
					err = nil
					break
				}
				return nil, errDecorate(err, "GroRead")
			*/
		}
		Coords = append(Coords, tmpcoords)
	}
	returned, err := NewMolecule(Coords, top, nil)
	//	fmt.Println("2 return!", top.Atom(1), returned.Coords[0].VecView(2)) ///////////////////////
	return returned, errDecorate(err, "GroRead")
}

func groReadSnap(gro *bufio.Reader, ReadTopol bool) (*v3.Matrix, []*Atom, error) {
	nm2A := 10.0
	chains := "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	line, err := gro.ReadString('\n') //we don't care about this line,but it has to be there
	if err != nil {
		return nil, nil, CError{fmt.Sprintf("Empty gro File: %s", err.Error()), []string{"bufio.Reader.ReadString", "groReadSnap"}}
	}
	line, err = gro.ReadString('\n')
	if err != nil {
		return nil, nil, CError{fmt.Sprintf("Malformed gro File: %s", err.Error()), []string{"bufio.Reader.ReadString", "groReadSnap"}}
	}

	natoms, err := strconv.Atoi(strings.TrimSpace(line))
	if err != nil {
		return nil, nil, CError{fmt.Sprintf("Wrong header for a gro file %s", err.Error()), []string{"strconv.Atoi", "groReadSnap"}}
	}
	var molecule []*Atom
	if ReadTopol {
		molecule = make([]*Atom, 0, natoms)
	}
	coords := make([]float64, 0, natoms*3)
	prevres := 0
	chainindex := 0
	for i := 0; i < natoms; i++ {
		line, err = gro.ReadString('\n')
		if err != nil {
			return nil, nil, CError{fmt.Sprintf("Failure to read gro File: %s", err.Error()), []string{"bufio.Reader.ReadString", "groReadSnap"}}
		}
		fields := strings.Fields(line)
		if len(fields) < 4 {
			break //meaning this line contains the unit cell vectors, and it is the last line of the snapshot
		}
		if ReadTopol {
			atom, c, err := read_gro_line(line)
			if err != nil {
				return nil, nil, CError{fmt.Sprintf("Failure to read gro File: %s", err.Error()), []string{"bufio.Reader.ReadString", "groReadSnap"}}
			}
			if atom.MolID < prevres {
				chainindex++
			}
			prevres = atom.MolID
			if chainindex >= len(chains) {
				chainindex = 0 //more chains inthe molecule than letters in the alphabet!
			}
			atom.Chain = string(chains[chainindex])
			molecule = append(molecule, atom)
			coords = append(coords, c...)
			//	fmt.Println(atom, c) //////////////////
			continue

		}
		c := make([]float64, 3, 3)
		for i := 0; i < 3; i++ {
			c[i], err = strconv.ParseFloat(strings.TrimSpace(line[20+(i*8):28+(i*8)]), 64)
			if err != nil {
				return nil, nil, err
			}
			c[i] = c[i] * nm2A //gro uses nm, goChem uses A.
		}
		coords = append(coords, c...)

	}
	mcoords, err := v3.NewMatrix(coords)
	//	fmt.Println(molecule) //, mcoords) ////////////////////////
	return mcoords, molecule, nil
}

//read_gro_line Parses a valid ATOM or HETATM line of a PDB file, returns an Atom
// object with the info except for the coordinates and b-factors, which  are returned
// separately as an array of 3 float64 and a float64, respectively
func read_gro_line(line string) (*Atom, []float64, error) {
	coords := make([]float64, 3, 3)
	atom := new(Atom)
	nm2A := 10.0
	var err error
	atom.MolID, err = strconv.Atoi(strings.TrimSpace(line[0:5]))
	if err != nil {
		return nil, nil, err
	}
	atom.MolName = strings.TrimSpace(line[5:10])
	atom.MolName1 = three2OneLetter[atom.MolName]
	atom.Name = strings.TrimSpace(line[10:15])
	atom.ID, err = strconv.Atoi(strings.TrimSpace(line[15:20]))
	//	fmt.Printf("%s|%s|%s|%s|\n", line[0:5], line[5:10], line[10:15], line[15:20]) ////////////
	if err != nil {
		return nil, nil, err
	}
	for i := 0; i < 3; i++ {
		coords[i], err = strconv.ParseFloat(strings.TrimSpace(line[20+(i*8):28+(i*8)]), 64)
		if err != nil {
			return nil, nil, err
		}
		coords[i] = coords[i] * nm2A //gro uses nm, goChem uses A.
	}

	atom.Symbol, _ = symbolFromName(atom.Name)
	return atom, coords, nil
}

//GoFileWrite writes the molecule described by mol and Coords into a file in the Gromacs
//gro format. If Coords has more than one elements, it will write a multi-state file.
func GroFileWrite(outname string, Coords []*v3.Matrix, mol Atomer) error {
	out, err := os.Create(outname)
	if err != nil {
		return CError{"Failed to write open file" + err.Error(), []string{"os.Create", "GroFileWrite"}}
	}
	defer out.Close()
	for _, v := range Coords {
		err := GroSnapWrite(v, mol, out)
		if err != nil {
			return errDecorate(err, "GoFileWrite")
		}
	}
	return nil
}

//GroSnapWrite writes a single snapshot of a molecule to an io.Writer, in the Gro format.
func GroSnapWrite(coords *v3.Matrix, mol Atomer, out io.Writer) error {
	A2nm := 0.1
	iowriterError := func(err error) error {
		return CError{"Failed to write in io.Writer" + err.Error(), []string{"io.Writer.Write", "GroSnapWrite"}}
	}
	if mol.Len() != coords.NVecs() {
		return CError{"Ref and Coords dont have the same number of atoms", []string{"GroSnapWrite"}}
	}
	c := make([]float64, 3, 3)
	_, err := out.Write([]byte(fmt.Sprintf("Written with goChem :-)\n%-4d\n", mol.Len())))
	if err != nil {
		return iowriterError(err)
	}
	//towrite := Coords.Arrays() //An array of array with the data in the matrix
	for i := 0; i < mol.Len(); i++ {
		//c := towrite[i] //coordinates for the corresponding atoms
		c = coords.Row(c, i)
		at := mol.Atom(i)
		//velocities are set to 0
		temp := fmt.Sprintf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", at.MolID, at.MolName, at.Name, at.ID, c[0]*A2nm, c[1]*A2nm, c[2]*A2nm, 0.0, 0.0, 0.0)
		_, err := out.Write([]byte(temp))
		if err != nil {
			return iowriterError(err)
		}

	}
	//the box vectors at the end of the snappshot
	_, err = out.Write([]byte("0.0 0.0 0.0\n"))
	if err != nil {
		return iowriterError(err)
	}
	return nil
}
