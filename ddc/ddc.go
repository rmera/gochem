/*
 * dcd.go, part of gochem
 *
 * Copyright 2012 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License  as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

/**********
The DDCMD trajectory format has 2 parts.

The first one, or "header" contains utf-8 text.
The header itself is surrounded by {} and contains a set of key=value pairs separated by ";"
some "value" fields can contain several "words" separated by spaces.
The header is parsed here into a "header" structure with one variable per "key", which
takes the corresponding "value" of that key. When the value consists of more than one
"word", the textis split by fields (i.e. words) and each word becomes an element of a
slice. The slice is assigned to the variable. The "keys" are all strings, but only some
of the "values" are meant to be strings or  slices of string. Some are to be parsed as
int or []int or as float (in goChem's case, float64).

The second part contains the coordinates and other data. Of these, goChem only recovers
coordinates. This second part may be in several formats, including text and binary.
goChem, as of now, only support binary trajectories.

***********/

package ddc

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
	"unicode"
	"unicode/utf8"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

type header struct {
	offset      int
	datatype    string
	l           [3]float64
	nfields     int
	nfiles      int
	nrecord     int
	fieldName   []string
	fieldType   []string
	time        float64 //Whatever is in the file must be multiplied by 0.001 for some unit reason
	species     []string
	fieldUnit   []string
	fieldFormat []string
}

func (H *header) parseObj(obj map[string]string) error {
	//This is tricky, and, maybe, not the best way to do it, but it will have to do, at least for now.
	//At least, the ugliness is hidden here.
	req := []string{"datatype", "h", "field_names", "field_types", "nfields", "nfiles", "nrecord", "time"}
	for _, v := range req {
		_, ok := obj[v]
		if !ok {
			return Error{"Oject lacks one or more required fields!: " + v, "", []string{"parseObj"}, true}
		}
	}
	H.datatype = obj["datatype"]
	hfields := strings.Fields(obj["h"])
	if len(hfields) < 9 {
		return Error{"field 'h' has less than the required 9 sub-fields", "", []string{"parseObj"}, true}
	}
	ls, err := parseNFloats(hfields[0], hfields[4], hfields[8])
	if err != nil {
		return errDecorate(err, "parseObj: Parsing h")
	}
	H.l = [3]float64{ls[0], ls[1], ls[2]}

	intf, err := parseNInts(obj["nfields"], obj["nfiles"], obj["nrecord"])
	if err != nil {
		return errDecorate(err, "parseObj: Parseing nfields, nfiles and nrecord")
	}
	H.nfields = intf[0]
	H.nfiles = intf[1]
	H.nrecord = intf[2]
	H.fieldName = strings.Fields(obj["field_names"])
	H.fieldType = strings.Fields(obj["field_types"])
	time, err := strconv.ParseFloat(strings.Fields(obj["time"])[0], 64)
	if err != nil {
		return Error{"Parsing time: " + err.Error(), "", []string{"parseObj"}, true}
	}
	const timeconv = 0.001 //yeah, I don't know what this is either.
	H.time = time * timeconv

	//Those were the required fields, now a few optionals
	if spec, ok := obj["species"]; ok {
		H.species = strings.Fields(spec)
	}
	if form, ok := obj["field_format"]; ok {
		H.fieldFormat = strings.Fields(form)
	}
	if units, ok := obj["field_units"]; ok {
		H.fieldUnit = strings.Fields(units)
	}

	return nil
}

//Container for an Charmm/NAMD binary trajectory file.
type DDCObj struct {
	h        *header
	binary   bool
	readable bool
	filename string
	ddc      *os.File
	lunit    float64 //lenght unit conversion factor
	natoms   int
	nfiles   int
	fieldFmt []string
	fieldBit []int
	endian   binary.ByteOrder
}

func (H *header) lenunits(D *DDCObj) error {
	//Maps the unit name with it's conversion factor to A (goChem's units for lenght)
	unitmap := map[string]float64{
		"A":    1,
		"Ang":  1,
		"Bohr": chem.Bohr2A,
		"a0":   chem.Bohr2A,
		"nm":   0.1,
		"um":   0.1 * 10E-3,
		"mm":   0.1 * 10E-6, //I wonder why are these units even supported.
	}
	found := false
	for i, v := range H.fieldName {
		if v != "rx" {
			continue
		}
		var ok bool
		unittext := H.fieldName[i]
		if D.lunit, ok = unitmap[unittext]; !ok {
			return Error{"Can't find unit: " + unittext, D.filename, []string{"lenunit"}, true}
		}
		found = true
		break
	}
	if !found {
		log.Printf("Couldn't find a section with lenght units in header. Will use the default (A)")
		D.lunit = 1.0
	}
	return nil

}

//New builds a new DCDObj object from a DCD trajectory file
func New(filename string) (*DDCObj, error) {
	traj := new(DDCObj)
	if err := traj.initRead(filename); err != nil {
		return nil, errDecorate(err, "New")
	}
	return traj, nil

}

//Readable returns true if the object is ready to be read from
//false otherwise. It doesnt guarantee that there is something
//to read.
func (D *DDCObj) Readable() bool {
	return D.readable
}

//This function parses a text string delimited by "{" and "}", which contains pairs formatted as key=value,
//where each pair is separated from the others by ";". The function returns a map of the keys and values,
//just as text, so it will be somebody else's job to parse the values to int, float64, etc.
func (D *DDCObj) parseObj(obj string) (map[string]string, error) {
	obj = strings.Replace(obj, "\n", "", -1)
	f := func(r rune) bool {
		if r == '{' || r == '}' {
			return true
		}
		return false
	}
	obj2 := strings.FieldsFunc(obj, f)[1]
	if obj2 == obj {
		return nil, Error{"Object not delimited by {}", D.filename, []string{}, true}
	}
	mapaq := map[string]string{} //So, I have to parse the items into numbers or whatever, later.
	listaq := strings.Split(obj2, ";")
	trim := strings.TrimSpace //just convenience
	for _, v := range listaq {
		kv := strings.Split(v, "=")
		mapaq[trim(kv[0])] = trim(kv[1])
	}

	return mapaq, nil
}

func (D *DDCObj) bitFormat() error {
	ft := D.h.fieldType
	fform := make([]string, 0, len(ft))
	fbits := make([]int, 0, len(ft))
	for _, v := range ft {
		bit, err := strconv.Atoi(string(v[1])) //This will not work with non-ascii values! which shouldn't happen, anyway.
		if err != nil {
			return Error{"Wrong bit format: " + v, D.filename, []string{"bitFormat"}, true}
		}
		if v[0] == 'u' && bit == 4 {
			fform = append(fform, "unit32") // L
		} else if v[0] == 'u' && bit == 8 {
			fform = append(fform, "uint64") //Q
		} else if v[0] == 'f' && bit == 4 {
			fform = append(fform, "float32") //f
		} else if v[0] == 'f' && bit == 8 {
			fform = append(fform, "float64") //d
		} else {
			return Error{"Bit format not supported: " + v, D.filename, []string{"bitFormat"}, true}
		}
		D.endian = binary.LittleEndian //They are all little endian!
		fbits = append(fbits, bit)
	}
	D.fieldFmt = fform
	D.fieldBit = fbits
	return nil

}

//initRead initializes a XtcObj for reading.
//It requires only the filename, which must be valid.
//It support big and little endianness, charmm or (namd>=2.1) and no
//fixed atoms.
func (D *DDCObj) initRead(name string) error {
	headerstr := ""
	offset := 0
	//The file starts with an utf-8 header, then continues with the binary part.
	//we'll read the header now.
	var doneheader bool

	headerf := bufio.NewReader(D.ddc)
	for {
		lineb, err := headerf.ReadBytes('\n')
		if err != nil {
			break
		}
		if !utf8.Valid(lineb) {
			break
		}
		line := string(lineb)
		if !isallspaces(line) && doneheader {
			break
		}
		if strings.HasPrefix(line, "}") {
			doneheader = true
		}
		headerstr += line
		offset += len([]byte(line)) //I am not even sure about this!
	}

	hstruct := new(header)
	hstruct.offset = offset
	obj, err := D.parseObj(headerstr)
	if err != nil {
		return errDecorate(err, "initRead: Trying to parse the header string")
	}
	err = hstruct.parseObj(obj)
	if err != nil {
		return errDecorate(err, "initRead: Trying to parse the header dictionary")
	}
	err = hstruct.parseObj(obj)
	if err != nil {
		return errDecorate(err, "initRead: Trying to parse the header dictionary")
	}
	err = hstruct.lenunits(D)
	if err != nil {
		return errDecorate(err, "initRead: Trying to find units")
	}
	if hstruct.datatype == "FIXRECORDBINARY" {
		D.binary = true
		err = D.bitFormat()
		if err != nil {
			return errDecorate(err, "initRead")
		}

	} else {
		D.binary = false
		//At least for now. I really don't have the energy to implement every possible sub-format.
		return Error{"goChem only support binary DDCMD trajectories", D.filename, []string{"initRead"}, true}

	}
	D.readable = true
	D.h = hstruct
	//We now check that all fields needed to read the trajectory are present
	tf := []string{"id", "rx", "ry", "rz"}
	for _, v := range tf {
		if !isInString(D.h.fieldName, v) {
			return Error{"File lacks at least one required term: " + v, D.filename, []string{"initRead"}, true}
		}
	}
	D.nfiles = D.h.nfiles
	return nil

}

func (D *DDCObj) nextbin(coord *v3.Matrix) {
	filenames := make([]string, 0, 3)
	basename := strings.Split(D.filename, "#")[0]
	for i := 0; i < D.nfiles; i++ {
		numbering := fmt.Sprintf("%06d", strconv.Itoa(i))
		newname := basename + "#" + numbering
		filenames = append(filenames, newname)
	}

}

//Next Reads the next frame in a DcDObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DDCObj) Next(keep *v3.Matrix) error {
	if !D.readable {
		return Error{TrajUnIni, D.filename, []string{"Next"}, true}
	}
	if D.dcdFields == nil {
		D.dcdFields = make([][]float32, 3, 3)
		D.dcdFields[0] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[1] = make([]float32, int(D.natoms), int(D.natoms))
		D.dcdFields[2] = make([]float32, int(D.natoms), int(D.natoms))
	}
	if err := D.nextRaw(D.dcdFields); err != nil {
		return errDecorate(err, "Next")
	}
	if keep == nil {
		return nil
	}
	if r, _ := keep.Dims(); int32(r) < D.natoms {
		panic("Not enough space in matrix")
	}
	//outBlock := make([]float64, int(D.natoms)*3, int(D.natoms)*3)
	for i := 0; i < int(D.natoms); i++ {
		k := i - i*(i/int(D.natoms))
		keep.Set(k, 0, float64(D.dcdFields[0][k]))
		keep.Set(k, 1, float64(D.dcdFields[1][k]))
		keep.Set(k, 2, float64(D.dcdFields[2][k]))
	}
	//	final := NewVecs(outBlock, int(D.natoms), 3)
	//	fmt.Print(final)/////////7
	return nil
}

//Next Reads the next frame in a XtcObj that has been initialized for read
//With initread. If keep is true, returns a pointer to matrix.DenseMatrix
//With the coordinates read, otherwiser, it discards the coordinates and
//returns nil.
func (D *DCDObj) nextRaw(blocks [][]float32) error {
	if len(blocks[0]) != int(D.natoms) || len(blocks[1]) != int(D.natoms) || len(blocks[2]) != int(D.natoms) {
		return Error{NotEnoughSpace, D.filename, []string{"nextRaw"}, true}
	}
	D.new = false
	if D.readLast {
		D.readable = false
		return newlastFrameError(D.filename, "nextRaw")
	}
	//if there is an extra block we just skip it.
	//Sadly, even when there is an extra block, it is not present in all
	//snapshots for some trajectories, so we must use the block size to see if
	//there is an extra block or if the X block starts inmediately
	var blocksize int32
	if D.extrablock {
		if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
			return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}
		}
		//If the blocksize is 4*natoms it means that the block is not an
		//extra block, but the X coordinates, and thus we must skip the following
		if blocksize != D.natoms*4 {
			if _, err := D.readByteBlock(blocksize); err != nil {
				return err
			}
			blocksize = 0
		}
	}
	//now get the coords, each as a slice of float32
	//X
	//we collect the X block size again only if it has not been collected before
	if blocksize == 0 {
		if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
			return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}

		}
	}
	err := D.readFloat32Block(blocksize, blocks[0])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//	fmt.Println("X", len(xblock)) //, xblock)
	//	fmt.Println("X", blocks[0])  ///////////////////////////////
	//Y
	//Collect the size first, then the rest
	if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
		return err
	}
	err = D.readFloat32Block(blocksize, blocks[1])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//	fmt.Println("Y", blocks[1])
	//Z
	if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}
	}
	err = D.readFloat32Block(blocksize, blocks[2])
	if err != nil {
		return errDecorate(err, "nextRaw")
	}
	//	fmt.Println("Z", blocks[2])
	//we skip the 4-D values if they exist. Apparently this is not present in the
	//last snapshot, so we use an EOF here to signal that we have read the last snapshot.
	if D.charmm && D.fourdim {
		if err := binary.Read(D.dcd, D.endian, &blocksize); err != nil {
			if err.Error() == "EOF" {
				D.readLast = true
				//	fmt.Println("LAST!")
			} else {
				return Error{err.Error(), D.filename, []string{"binary.Read", "nextRaw"}, true}
			}
		}
		if !D.readLast {
			if _, err := D.readByteBlock(blocksize); err != nil {
				return errDecorate(err, "nextRaw")
			}
		}
	}
	return nil

}

//Queries the size of a block, and reads its contents into block, which must have the
//appropiate size.
func (D *DCDObj) readFloat32Block(blocksize int32, block []float32) error {
	var check int32
	//	fmt.Println("blockf",blocksize)
	if err := binary.Read(D.dcd, D.endian, block); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Read", "readFloat32Block"}, true}

	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return Error{err.Error(), D.filename, []string{"binary.Read", "readFloat32Block"}, true}
	}
	if check != blocksize {
		return Error{WrongFormat, D.filename, []string{"readFloat32Block"}, true}
	}
	return nil
}

//Queries the size of a block, make a slice of a quarter of that size
//and reads that amount of float32. This function is used for the
//
func (D *DCDObj) readByteBlock(blocksize int32) ([]byte, error) {
	var check int32
	//	fmt.Println("blockb",blocksize)
	block := make([]byte, blocksize, blocksize)
	if err := binary.Read(D.dcd, D.endian, block); err != nil {
		return nil, Error{err.Error(), D.filename, []string{"binary.Read", "readByteBlock"}, true}

	}
	if err := binary.Read(D.dcd, D.endian, &check); err != nil {
		return nil, Error{err.Error(), D.filename, []string{"binary.Read", "readByteBlock"}, true}

	}
	if check != blocksize {
		return nil, Error{SecurityCheckFailed, D.filename, []string{"readByteBlock"}, true}
	}
	return block, nil
}

//Len returns the number of atoms per frame in the XtcObj.
//XtcObj must be initialized. 0 means an uninitialized object.
func (D *DDCObj) Len() int {
	return int(D.natoms)
}

//Helper functions

//adds as many "0" as needed for an ASCII string to be of length target
func zfill(target int, tofill string) string {
	l := len([]byte(tofill))
	for ; l <= target; l++ {
		tofill = "0" + tofill //yeah, not efficient. I suppose "target" will never be greater than ~10 so this is fine.
		//if we find out it's not, well, changing the internal working of this function is easy.
	}
	return tofill

}

func isallspaces(test string) bool {
	testr := []rune(test)
	for _, v := range testr {
		if !unicode.IsSpace(v) {
			return false
		}
	}
	return true
}

func parseNFloats(s ...string) ([]float64, error) {
	ret := make([]float64, 0, len(s))
	for _, v := range s {
		f, err := strconv.ParseFloat(v, 64)
		if err != nil {
			return nil, Error{"Couldn't parse a float", "", []string{"parseNFloats"}, true}
		}
		ret = append(ret, f)
	}
	return ret, nil
}
func parseNInts(s ...string) ([]int, error) {
	ret := make([]int, 0, len(s))
	for _, v := range s {
		f, err := strconv.Atoi(v)
		if err != nil {
			return nil, Error{"Couldn't parse an int", "", []string{"parseNints"}, true}
		}
		ret = append(ret, f)
	}
	return ret, nil
}

//I guess this will all go after Go 1.18
func isInString(container []string, test string) bool {
	if container == nil {
		return false
	}
	for _, i := range container {
		if test == i {
			return true
		}
	}
	return false
}

//Errors

//errDecorate is a helper function that asserts that the error is
//implements chem.Error and decorates the error with the caller's name before returning it.
//if used with a non-chem.Error error, it will just return the error..
func errDecorate(err error, caller string) error {
	err2, ok := err.(chem.Error) //I know that is the type returned byt initRead
	if !ok {
		return err
	}
	err2.Decorate(caller)
	return err2
}

//Error is the general structure for DCD trajectory errors. It fullfills  chem.Error and chem.TrajError
type Error struct {
	message  string
	filename string //the input file that has problems, or empty string if none.
	deco     []string
	critical bool
}

func (err Error) Error() string {
	return fmt.Sprintf("dcd file %s error: %s", err.filename, err.message)
}

//Decorate Adds new information to the error
func (E Error) Decorate(deco string) []string {
	//Even thought this method does not use a pointer as a receiver, and tries to alter the received,
	//it should work, since E.deco is a slice, and hence a pointer itself.

	if deco != "" {
		E.deco = append(E.deco, deco)
	}
	return E.deco
}

//Filename returns the file to which the failing trajectory was associated
func (err Error) FileName() string { return err.filename }

//Format returns the format of the file (always "dcd") associated to the error
func (err Error) Format() string { return "dcd" }

//Critical returns true if the error is critical, false otherwise
func (err Error) Critical() bool { return err.critical }

const (
	TrajUnIni           = "Traj object uninitialized to read"
	ReadError           = "Error reading frame"
	UnableToOpen        = "Unable to open file"
	SecurityCheckFailed = "FailedSecurityCheck"
	WrongFormat         = "Wrong format in the DCD file or frame"
	NotEnoughSpace      = "Not enough space in passed blocks"
	EOF                 = "EOF"
)

//lastFrameError implements chem.LastFrameError
type lastFrameError struct {
	deco     []string
	fileName string
}

//lastFrameError does nothing
func (E lastFrameError) NormalLastFrameTermination() {}

func (E lastFrameError) FileName() string { return E.fileName }

func (E lastFrameError) Error() string { return "EOF" }

func (E lastFrameError) Critical() bool { return false }

func (E lastFrameError) Format() string { return "dcd" }

func (E lastFrameError) Decorate(deco string) []string {
	//Even thought this method does not use a pointer as a receiver, and tries to alter the received,
	//it should work, since E.deco is a slice, and hence a pointer itself.
	if deco != "" {
		E.deco = append(E.deco, deco)
	}
	return E.deco
}

func newlastFrameError(filename string, caller string) *lastFrameError {
	e := new(lastFrameError)
	e.fileName = filename
	e.deco = []string{caller}
	return e
}
