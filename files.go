/*
 * files.go, part of gochem.
 * 
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 * 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
import  "github.com/skelterjohn/go.matrix"


//Pdb_read family

	

//A map for assigning mass to elements. 
//Note that just common "bio-elements" are present
var symbolMass = map[string] float64 {
    "H":  1.0,
    "C":  12.01,
    "O": 16.00,
    "N": 14.01,
    "P": 30.97,
    "S": 32.06,
    "Se": 78.96,
    "K": 39.1,
    "Ca": 40.08,
    "Mg": 24.30,
    "Cl": 35.45,
    "Na": 22.99,
    "Cu":  63.55,
    "Zn":  65.38,
    "Co": 58.93,
    "Fe": 55.84,
    "Mn": 54.94,   
	}

//A map between 3-letters name for aminoacidic residues to the corresponding 1-letter names.
var three2OneLetter = map[string] byte{
	"SER": 'S',
	"THR": 'T',
	"ASN": 'N',
	"GLN": 'Q',
	"SEC": 'U',  //Selenocysteine!
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
func symbolFromName(name string) (string, error){
	symbol:=""
	if len(name)==4 || name[0]=='H'{  //I thiiink only Hs can have 4-char names in amber.
		symbol="H"
		}else if name[0]=='C'{ //Ca is not considered here
		if name=="CU"{
			symbol="Cu"
			}else if name=="CO"{
			symbol="Co"
			}else if name=="CL"{
			symbol="Cl"
			}else{
			symbol="C"
			}
		}else if name[0]=='N'{
		if name=="NA"{
			symbol="Na"
			}else{
			symbol="N"
			}
		}else if name[0]=='O'{
		symbol="O"
		}else if name[0]=='P'{
		symbol="P"	
		}else if name[0]=='S'{
		if name=="SE"{
			symbol="Se"
			}else{
			symbol="S"
			}
		}else if name[0:2]=="ZN"{
			symbol="Zn"
			}
	if symbol==""{
		return symbol, fmt.Errorf("Couldn't guess symbol from PDB name")
		}
	return symbol, nil 
	}


//Parses a valid ATOM or HETATM line of a PDB file, returns an Atom
// object with the info except for the coordinates and b-factors, which  are returned
// separately as an array of 3 float64 and a float64, respectively
func read_full_pdb_line(line string, read_additional bool, contlines int) (*Atom, []float64, float64,error) {
	err:=make([]error,7,7) //accumulate errors to check at the end of the readed line.
	coords:=make([]float64,3,3)
	atom:=new(Atom)
	atom.Het=strings.HasPrefix(line, "HETATM")  //this is called twice in the worst case. should fix
	atom.Id,err[0]=strconv.Atoi(strings.TrimSpace(line[6:12]))
	atom.Name=strings.TrimSpace(line[12:16])
	//PDB says that pos. 17 is for other thing but I see that is 
	//used for residue name in many cases*/
	atom.Molname=line[17:20]
	atom.Molname1=three2OneLetter[atom.Molname]
	atom.Chain=line[21] //currently this is read to a byte
	atom.Molid,err[1]=strconv.Atoi(strings.TrimSpace(line[22:26]))
	//Here we shouldn't need TrimSpace, but I keep it just in case someone
	 // doesn's use all the fields when writting a PDB*/
	coords[0],err[2]=strconv.ParseFloat(strings.TrimSpace(line[30:38]),64)
	coords[1],err[3]=strconv.ParseFloat(strings.TrimSpace(line[38:46]),64)
	coords[2],err[4]=strconv.ParseFloat(strings.TrimSpace(line[46:54]),64)
	atom.Occupancy,err[5]=strconv.ParseFloat(strings.TrimSpace(line[54:60]),64)
	//If I try not to declare this and just use :=, I get an "expected identifier" error
	var bfactor float64
	bfactor,err[6]=strconv.ParseFloat(strings.TrimSpace(line[60:66]),64)
	//we try to read the additional only if indicated and if it is there
	// In this part we don't catch errors. If something is missing we 
	// just ommit it
	if read_additional && len(line)>=80{
		atom.Symbol=(strings.TrimSpace(line[76:78]))
		atom.Charge=float64(line[78]) //strconv.ParseFloat(strings.TrimSpace(line[78:78]),64)
		if strings.Contains(line[79:79],"-"){
			atom.Charge=-1.0*atom.Charge 
			}
		}
		
	//This part tries to guess the symbol from the atom name, if it has not been read 
	//No error checking here, just fills symbol with the empty string the function returns
	if len(atom.Symbol)==0{
		atom.Symbol,_=symbolFromName(atom.Name)
		}
		
		
	for i:=range(err){
		if err[i]!=nil{
			//Here I should add the line number to the returned error.
			return nil, nil, 0, err[i]
			}
		}
	if atom.Symbol!=""{
		atom.Mass=symbolMass[atom.Symbol] //Not error checking
		}
	return atom, coords, bfactor, nil
	}
	
/*Parses a PDB line if only the coordinates and bfactors are to be read. right now the returned error
 * is always nil the funcion just cancels the execution of the program if there are trouble
 * This might change in the future */	
func read_onlycoords_pdb_line(line string, contlines int) ([]float64,float64, error) {
	coords:=make([]float64,3,3)
	err:=make([]error,4,4)
	var bfactor float64 //I dont get why I must declare this instead of using :=
	//I get an "expected identifier" error if I do so.
	coords[0],err[0]=strconv.ParseFloat(strings.TrimSpace(line[30:38]),64)
	coords[1],err[1]=strconv.ParseFloat(strings.TrimSpace(line[38:46]),64)
	coords[2],err[2]=strconv.ParseFloat(strings.TrimSpace(line[46:54]),64)
	bfactor, err[3] = strconv.ParseFloat(strings.TrimSpace(line[60:66]),64)
	//this will take care of any error
	for i:=range(err){
		if err[i]!=nil{
			//Here I should add the line number to the returned error.
			return nil, 0, err[i]
			}
		}
	return coords,bfactor,nil
	}


//PdbRead reads the atomic entries for a PDB file, returns a bunch of without coordinates,
// and the coordinates in a separate array of arrays. If there is one frame in the PDB
// the coordinates array will be of lenght 1. It also returns an error which is not 
// really well set up right now.
func PdbRead(pdbname string, read_additional bool) ([]*Atom, []*matrix.DenseMatrix,[][]float64, error){
	molecule:=make([]*Atom,0) //I thiiink is more efficient to have pointers here
	coords:=make([][]float64,1,1)
	coords[0]=make([]float64,0)
	bfactors:=make([][]float64,1,1)
	bfactors[0]=make([]float64,0)
	first_model:=true //are we reading the first model? if not we only save coordinates
	pdbfile, err := os.Open(pdbname)
	if err!= nil{
		//fmt.Println("Unable to open file!!")
		return molecule,nil,nil,err 
		}
	defer pdbfile.Close()
	pdb := bufio.NewReader(pdbfile)
	contlines:=1 //count the lines read to better report errors
	for {
		line, err := pdb.ReadString('\n')
		if err != nil {
			fmt.Println("PDB reading complete")  /***change this to stderr************/
			break
		contlines++ //count all the lines even if empty.
			}
		if len(line)<4{
			continue
			}
		//here we start actually reading
		/**There might be a bug for not copying the string (symbol, name, etc) but just assigning the slice
		 * which is a reference. Check!**/
		var c=make([]float64,3,3)
		var bfactemp float64 //temporary bfactor
		var atomtmp *Atom
	//	var foo string // not really needed
		if strings.HasPrefix(line, "ATOM") || strings.HasPrefix(line, "HETATM"){
			if !first_model{
				c,bfactemp,err=read_onlycoords_pdb_line(line, contlines)
				}else{
				atomtmp=new(Atom)
				atomtmp,c,bfactemp,err=read_full_pdb_line(line, read_additional, contlines)
				if err!=nil{
					return molecule,nil,nil,err 
					}
				//atom data other than coords is the same in all models so just read for the first.
				molecule=append(molecule,atomtmp)
				}
			//coords are appended for all the models
			//we add the coords to the latest frame of coordinaates
			coords[len(coords)-1]=append(coords[len(coords)-1],c[0],c[1],c[2])
			bfactors[len(bfactors)-1]=append(bfactors[len(bfactors)-1],bfactemp)
			} else if strings.HasPrefix(line,"MODEL"){
			modelnumber:=1  //apparently in PDBs the counts starts from 1 
			modelnumber,_=strconv.Atoi(strings.TrimSpace(line[6:]))
			if modelnumber > 1{
				first_model=false
				coords=append(coords,make([]float64,0)) //new bunch of coords for a new frame
				bfactors=append(bfactors,make([]float64,0))
				}
			}
		}	
	//This could be done faster if done in the same loop where the coords are read
	//Instead of having another loop just for them.
	frames:=len(coords)
	mcoords:=make([]*matrix.DenseMatrix,frames,frames) //Final thing to return
	for i:=0;i<frames;i++{
		mcoords[i]=matrix.MakeDenseMatrix(coords[i],len(coords[i])/3,3)
		}
	/**tests for debugging**/
	fmt.Println("Coords read", mcoords[0].GetMatrix(0,0,3,3),"Atoms: ",len(molecule)," Coords: ",mcoords[0].Rows()) 
	//We ensure an slice (not a copy!) is passed. should save memory and cpu
	return molecule[:], mcoords[:],bfactors[:], err
	}
//End Pdb_read family


//PdbWrite writes a PDB file for the molecule mol with the file name pdbname.
func PdbWrite(mol *Molecule, pdbname string) error{
	err:=mol.Corrupted()
	if err!=nil{
		return err
		}
	out, err := os.Create(pdbname)
	if err!=nil{
		return err
		}
	defer out.Close()
	fmt.Fprint(out,"REMARK     WRITTEN WITH GOCHEM :-)\n")
	for j:= range mol.Coords{
		towrite:=mol.Coords[j].Arrays()  //An array of array with the data in the matrix
		chainprev:=mol.Atoms[0].Chain  //this is to know when the chain changes.
		fmt.Fprintf(out, "MODEL %d\n",j)
//		fmt.Println("chain", mol.Atoms[0])		
		for i:=range mol.Atoms{
			if mol.Atoms[i].Chain!=chainprev{
				fmt.Fprintln(out,"TER")
				chainprev=mol.Atoms[i].Chain
				}
			first:="ATOM"
			if mol.Atoms[i].Het{
				first="HETATM"
				}
			//Corrupted will check for broken bfactors and complete with zeroes
			//So no need to check here.
			c:=towrite[i] //coordinates for the corresponding atoms
			if len(mol.Atoms[i].Name)<4{
		//		fmt.Println("chain", mol.Atoms[i])
				_,err=fmt.Fprintf(out,"%-6s%5d  %-3s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",first, mol.Atoms[i].Id, mol.Atoms[i].Name, mol.Atoms[i].Molname, mol.Atoms[i].Chain,
							mol.Atoms[i].Molid, c[0],c[1],c[2], mol.Atoms[i].Occupancy, mol.Bfactors[j][i], mol.Atoms[i].Symbol)
				//4 chars for the atom name are used when hydrogens are included.	
				//This has not been tested
				}else if len(mol.Atoms[i].Name)==4{
				_,err=fmt.Fprintf(out,"%-6s%5d %4s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",first, mol.Atoms[i].Id, mol.Atoms[i].Name, mol.Atoms[i].Molname, mol.Atoms[i].Chain,
							mol.Atoms[i].Molid, c[0],c[1],c[2], mol.Atoms[i].Occupancy, mol.Bfactors[j][i], mol.Atoms[i].Symbol)			
				}else{
				err=fmt.Errorf("Cant print PDB line")
				}
			if err!=nil{
				return err
				}
			}
		fmt.Fprint(out,"ENDMDL\n")
		}
	fmt.Fprint(out,"END\n")
	return nil
	}
	
	

//XyzRead reads an xyz file, returns a slice of Atom objects, and slice of matrix.DenseMatrix and an error.
func XyzRead(xyzname string,) ([]*Atom, []*matrix.DenseMatrix, error){
	xyzfile, err:= os.Open(xyzname)
	if err!= nil{
		//fmt.Println("Unable to open file!!")
		return nil,nil,err 
		}
	defer xyzfile.Close()
	xyz := bufio.NewReader(xyzfile)	
	line, err := xyz.ReadString('\n')
//	fmt.Println("line: ", line) /////////////////////////////////7
	if err != nil {
		return nil, nil, fmt.Errorf("Ill formatted XYZ file!")
		}
//	var natoms int
	natoms,err:=strconv.Atoi(strings.TrimSpace(line))
	if err != nil {
		return nil, nil, fmt.Errorf("Ill formatted XYZ file!")
		}
//	fmt.Println("natoms: ", natoms)///////////////////////////////7
	molecule:=make([]*Atom,natoms,natoms)
	coords:=make([]float64,natoms*3,natoms*3)
	//THIS COULD BECOME A BUG IF TRY TO READ AN EMPTY FILE!!!!!!!!!!!!!!!!!!!
	//MUST FIX!!!!!!!
	_,_=xyz.ReadString('\n') //We dont care about this line
//	var i int
//	var line string
	errs:=make([]error,3,3)
	for i:=0;i<natoms;i++{
		line, errs[0] = xyz.ReadString('\n')
		if errs[0] != nil {  //inefficient, (errs[1] can be checked once before), but clearer.
			break
			}
		fields:=strings.Fields(line)
		if len(fields)<4{
			errs[0]=fmt.Errorf("Line number %d in file %s ill formed",i,xyzname)
			break
			}
		molecule[i]=new(Atom)
		molecule[i].Symbol=fields[0]
		molecule[i].Mass=symbolMass[molecule[i].Symbol]
		coords[i*3],errs[0]=strconv.ParseFloat(fields[1],64)
		coords[i*3+1],errs[1]=strconv.ParseFloat(fields[2],64)
		coords[i*3+2],errs[2]=strconv.ParseFloat(fields[3],64)
	}	
	//This could be done faster if done in the same loop where the coords are read
	//Instead of having another loop just for them.
	for _,i:=range errs{
		if i!=nil{
			return nil, nil, i
			}
		}
	mcoords:=make([]*matrix.DenseMatrix,1,1)
	mcoords[0]=matrix.MakeDenseMatrix(coords,natoms,3)
	return molecule[:], mcoords[:], errs[0]
	}




//XyzWrite writes the frame frame of molecule mol in an XYZ file with name xyzname which will
//be created fot that. If the file exist it will be overwriten.
func XyzWrite(mol *Molecule, frame int, xyzname string) error{
	err:=mol.Corrupted()
	if err!=nil{
		return err
		}
	out, err := os.Create(xyzname)
	if err!=nil{
		return err
		}
	defer out.Close()
	fmt.Fprintf(out,"%-4d\n",len(mol.Atoms))
	fmt.Fprintf(out,"\n")
	towrite:=mol.Coords[frame].Arrays()  //An array of array with the data in the matrix	
	for i:=range mol.Atoms{
		c:=towrite[i] //coordinates for the corresponding atoms
		_,err=fmt.Fprintf(out,"%-2s  %8.3f%8.3f%8.3f \n",mol.Atoms[i].Symbol, c[0],c[1],c[2])
		if err!=nil{
			return err
			}
		}
	return nil
	}
