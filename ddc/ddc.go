/*
 * ddc.go, part of gochem
 *
 *
 * Copyright (c) 2021 Raul Mera Adasme <raul_dot_mera_changeforat_usach_dot_cl>
 *
 * Partially based on the implementation at:
 * https://github.com/XiaohuaZhangLLNL/mdanalysis
 * Which is (c) of the ddcMD MDAnalysis fork authors.
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
This reads ddcMD trajectory files

The ddcMD trajectory format has 2 parts.

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
The coordinates are in sets of numbers. All sets in a trajectory have the same
number/type of elements.
How many, what are, and what "format" (i.e. float64, uint32, etc.)
each element in a set is, is given by the fieldXXXX data from
the header. We should have at least these fields in a set: id, rx, ry and rz. We only
care about the 3 last ones, and discard id, and any other additional field.

We read sets one after the other, until the file is over.

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
	//	offset      int
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

//parseObj takes a map parsed from a ddcMD formatted header, and casts the values from strings to
//whatever type is needed, putting it into the receiver's respective fields.
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

//Container for a ddcMD trajectory file.
type DDCObj struct {
	h            *header
	offset       int64
	binary       bool
	readable     bool
	filename     []string //The filename for each frame (not working right now, and could change, but that's my plan right now).
	framecounter int
	ddc          *os.File // I think I should set this to point to the file and point to wherever we should start reading next.It would not be attached to one particular file as in my other trajectory formats.
	lunit        float64  //lenght unit conversion factor
	natoms       int
	nfiles       []int //how many files on each frame
	fieldFmt     []string
	fieldBit     []int
	endian       binary.ByteOrder
}

//Obtains the len unit conversion from the units in the DDC object to A
func (H *header) lenunits(D *DDCObj) error {
	//Maps the unit name with it's conversion factor to A (goChem's units for lenght)
	unitmap := map[string]float64{
		"A":    1,
		"Ang":  1,
		"Bohr": chem.Bohr2A,
		"a0":   chem.Bohr2A,
		"nm":   0.1,
		"um":   0.1 * 10E-3,
		"mm":   0.1 * 10E-6,
	}
	found := false
	for i, v := range H.fieldName {
		if v != "rx" {
			continue
		}
		var ok bool
		unittext := H.fieldUnit[i]
		if D.lunit, ok = unitmap[unittext]; !ok {
			return Error{"Can't find unit: " + unittext, D.filename[D.framecounter], []string{"lenunit"}, true}
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
		return nil, Error{"Object not delimited by {}", D.filename[D.framecounter], []string{}, true}
	}
	mapaq := map[string]string{} //I still have to parse the values into numbers or whatever.
	listaq := strings.Split(obj2, ";")
	trim := strings.TrimSpace //just convenience
	for _, v := range listaq {
		if isallspaces(v) {
			continue //it should be the same as "break" really, we should be done here.
		}
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
		bit, err := strconv.Atoi(string(v[1]))
		if err != nil {
			return Error{"Wrong bit format: " + v, D.filename[D.framecounter], []string{"bitFormat"}, true}
		}
		if v[0] == 'u' && bit == 4 {
			fform = append(fform, "uint32") // L
		} else if v[0] == 'u' && bit == 8 {
			fform = append(fform, "uint64") //Q
		} else if v[0] == 'f' && bit == 4 {
			fform = append(fform, "float32") //f
		} else if v[0] == 'f' && bit == 8 {
			fform = append(fform, "float64") //d
		} else {
			return Error{"Bit format not supported: " + v, D.filename[D.framecounter], []string{"bitFormat"}, true}
		}
		D.endian = binary.LittleEndian //They are all little endian
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
	var offset int64
	D.filename = append(D.filename, name)
	//The file starts with an utf-8 header, then continues with the binary part.
	//we'll read the header now.
	var doneheader bool
	ddc, err := os.Open(name)
	if err != nil {
		return Error{"Couldn't open first frame/header file ", name, []string{"initRead"}, true}
	}
	headerf := bufio.NewReader(ddc)
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
		offset += int64(len([]byte(line))) //This is for file.Seek() use later, so we skip all the header.
	}

	hstruct := new(header)
	D.offset = offset
	obj, err := D.parseObj(headerstr)
	if err != nil {
		return errDecorate(err, "initRead: Trying to parse the header string")
	}
	err = hstruct.parseObj(obj)
	if err != nil {
		return errDecorate(err, "initRead: Trying to parse the header dictionary")
	}
	D.h = hstruct
	err = D.h.lenunits(D)
	if err != nil {
		return errDecorate(err, "initRead: Trying to find units")
	}
	if D.h.datatype == "FIXRECORDBINARY" {
		D.binary = true
		err = D.bitFormat()
		if err != nil {
			return errDecorate(err, "initRead")
		}

	} else {
		D.binary = false
		//At least for now.
		return Error{"goChem only support binary DDCMD trajectories", D.filename[D.framecounter], []string{"initRead"}, true}

	}
	D.readable = true
	//We now check that all fields needed to read the trajectory are present
	//I actually don't need the "id" field. For now, I'll stick to the reference implementation and require it.
	tf := []string{"id", "rx", "ry", "rz"}
	for _, v := range tf {
		if !isInString(v, D.h.fieldName) {
			return Error{"File lacks at least one required term: " + v, D.filename[D.framecounter], []string{"initRead"}, true}
		}
	}
	D.nfiles = append(D.nfiles, D.h.nfiles)
	D.natoms = D.h.nrecord
	return nil

}

func (D *DDCObj) nextbin(coord *v3.Matrix) error {
	if coord == nil {
		//This will only work if each frame is in its own file!
		//otherwise I have to actually read the frame, and ignore it.
		//well, I could calculate the offset with the data
		//present per atom, and the number of atoms, all info I have.
		//Just adding all the fieldbit (which are actualy bytes) and
		//multiplying for the atoms plus the previous offset should give me the
		//right value to use f.Seek() later.
		//
		D.framecounter++
		return nil
	}
	if D.framecounter > 0 {
		D.initRead(D.filename[D.framecounter]) //All I actually need is to generate the offset value!
		//I could do something more efficient, i.e., that doesn't read all the header info again
	}
	filenames := make([]string, 0, 3)
	basename := strings.Split(D.filename[D.framecounter], "#")[0]
	//I should change this for a system where I just look for all files
	//named basename#SOMENUMBER and add it, ordered by the number.
	//I am kinda assuming that snapshots other than the first don't have a header,
	//because that would be crazy.
	for i := 0; i < D.nfiles[D.framecounter]; i++ {
		numbering := fmt.Sprintf("%06d", i)
		newname := basename + "#" + numbering
		filenames = append(filenames, newname)
	}
	atomsread := 0
	for i, file := range filenames {
		f, err := os.Open(file)
		if err != nil {
			return Error{fmt.Sprintf("Couldn't open frame's %d file", i), file, []string{"nextbin"}, true}

		}
		if i == 0 {
			f.Seek(D.offset, 0) //skip the header
		}
	ATOM: //This is the loop that reads coordinates for every atom, the loop inside it reads the x, y and z
		for err == nil {
			if atomsread >= D.natoms {
				break
			}
			coordtmp := [3]float64{}
			c := 0
			for i, format := range D.fieldFmt {
				var val float64
				name := D.h.fieldName[i]
				if !isInString(name, []string{"rx", "ry", "rz"}) {
					err = readrejected(f, D.endian, format)
					if err != nil {
						break ATOM
					}
				} else {
					val, err = readfloat(f, D.endian, format)
					if err != nil {
						break ATOM
					}
					coordtmp[c] = val * D.lunit
					c++ //ok this is kinda funny
				}
			}
			coord.Set(atomsread, 0, coordtmp[0])
			coord.Set(atomsread, 1, coordtmp[1])
			coord.Set(atomsread, 2, coordtmp[2])
			atomsread++

		}
		if err != nil {
			if atomsread == D.natoms && err.Error() == "EOF" {
				return newlastFrameError(file, "nextbin")
			} else {
				return Error{fmt.Sprintf("Unexpected error while reading coordinates: %s", err.Error()), file, []string{"nextbin"}, true}
			}
		}

	}

	D.offset = 0
	D.framecounter++
	return nil
}

//Next Reads the next frame in a DDcObj that has been initialized for read
//With initread. If keep is a v3.Matrix pointer, it will fill it with the values
//for the coordinates in the frame. If keep is nil, it will discard the frame.
func (D *DDCObj) Next(keep *v3.Matrix) error {
	if !D.readable {
		return Error{"Trajectory not ready to read", "", []string{"Next"}, true}
	}
	err := D.nextbin(keep)
	if err != nil && err.Error() != "EOF" {
		errDecorate(err, "Next")
	}
	return err
}

//Len returns the number of atoms per frame in the XtcObj.
//XtcObj must be initialized. 0 means an uninitialized object.
func (D *DDCObj) Len() int {
	return int(D.natoms)
}

//Helper functions

func readrejected(f *os.File, endian binary.ByteOrder, format string) error {
	var err error
	if strings.Contains(format, "uint") {
		err = readint(f, endian, format) //the value itself is discarded
	} else {
		_, err = readfloat(f, endian, format) //the value itself is discarded

	}
	return err

}

func readint(f *os.File, endian binary.ByteOrder, format string) error {
	var ret uint64
	var t uint32
	var err error
	if format == "uint32" {
		err = binary.Read(f, endian, &t)
	} else if format == "uint64" {
		err = binary.Read(f, endian, &ret)
	} else {
		return Error{fmt.Sprintf("Wrong format, should be uint64 or uint32, is %s", format), "", []string{"readint"}, true}
	}
	return err
}

func readfloat(f *os.File, endian binary.ByteOrder, format string) (float64, error) {
	var ret float64
	var t float32
	var err error
	if format == "float32" {
		err = binary.Read(f, endian, &t)
		ret = float64(t)
	} else if format == "float64" {
		err = binary.Read(f, endian, &ret)
	} else {
		return 0, Error{fmt.Sprintf("Wrong format, should be float64 or float32, is %s", format), "", []string{"readfloat"}, true}
	}
	return ret, err
}

func index(test string, where []string) int {
	for i, v := range where {
		if test == v {
			return i
		}
	}
	return -1
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
			return nil, Error{"Couldn't parse a float: " + v, "", []string{"parseNFloats"}, true}
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
			return nil, Error{"Couldn't parse an int: " + v, "", []string{"parseNints"}, true}
		}
		ret = append(ret, f)
	}
	return ret, nil
}

//I guess this will all go after Go 1.18
func isInString(test string, container []string) bool {
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
func (err Error) Format() string { return "ddc" }

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

func (E lastFrameError) Format() string { return "ddc" }

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
