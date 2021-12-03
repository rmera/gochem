package csf

import (
	"bufio"
	"bytes"
	"compress/flate"
	"compress/gzip"
	"compress/lzw"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"unsafe"

	//	"github.com/rmera/gochem"
	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

/***
Format spec:

The compatible and simple format is nothing but a compressed xyz

An STF file might have 3 file extension, which indicate
how they are compressed:
	stf: Compressed with deflate, at any level
	stz: Compressed with gzip (no idea why this is different than stf, but gunzip won't open my deflate-compressed files
	stl: Compressed with lzw

This spec will not establish what to do in case of a filename not ending with f, z or l, but the reference implementation
attempts to open anyway, and assumes it to be compressed with deflate. I repeat, that behavior is _not_ part of the spec, and might change.

A STF file may only contain ASCII symbols.

A STF file has a "header" starting in the first line, until a line that contains "**".
Each line of the header must be a pair key=value

After that, the file has one line per atom with 3 numbers (x y z coordinates), per frame
each frame ends with a line containing only "*".
**/

//Another option is make the trajectory a TAR file containing compressed snapshots. It would
//use somewhat more space, but I think it would be cheaper memory-wise to write.

const (
	lzwLitwidth int = 8
)

//Write!
type CsfW struct {
	f         *os.File
	h         io.WriteCloser
	natoms    int
	filename  string
	writeable bool
	symbols   []string
	header    map[string]string
	frames    int
}

func (S *CsfW) Close() {
	S.h.Close()
	S.f.Close()
	S.writeable = false
	return
}

func (S *CsfW) Len() int {
	return S.natoms
}

//compatibility with Gonum
func (S *CsfW) WNextDense(dcoord *mat.Dense) error {
	coord := v3.Dense2Matrix(dcoord)
	err := S.WNext(coord)
	if err != nil {
		err = errDecorate(err, "WNextDense")
	}
	return err
}

func (S *CsfW) WNext(coord *v3.Matrix) error {
	if !S.writeable {
		return Error{TrajUnIniWrite, S.filename, []string{"WNext"}, true}
	}
	if coord == nil {
		return Error{NilCoordinates, S.filename, []string{"WNext"}, true}
	}
	v := coord.NVecs()
	if v != S.natoms {
		return Error{fmt.Sprintf("%d coordinates given, but %d expected", v, S.natoms), S.filename, []string{"WNext"}, true}
	}

	S.h.Write([]byte(strconv.Itoa(S.natoms) + "\n"))
	//we write a header only in the first frame.
	if S.frames == 0 {
		str := make([]string, 0, len(S.header)+1) //one  more element is for the final "\n"
		for k, w := range S.header {
			str = append(str, fmt.Sprintf("%s=%s", k, w))
		}
		str = append(str, "\n")
		S.h.Write([]byte(strings.Join(str, "**")))
	} else {
		S.h.Write([]byte("\n")) //frames other than the first don't have anything here.
	}
	//now the coordinates
	strs := make([]string, 5)
	for i := 0; i < v; i++ {
		ff := strconv.FormatFloat
		prec := 3
		strs[0] = S.symbols[i]
		strs[1] = ff(coord.At(i, 0), 'f', prec, 64)
		strs[2] = ff(coord.At(i, 1), 'f', prec, 64)
		strs[3] = ff(coord.At(i, 2), 'f', prec, 64)
		strs[4] = "\n"
		str := strings.Join(strs, " ")
		//str := fmt.Sprintf("%07.3f %07.3f %07.3f\n", coord.At(i, 0), coord.At(i, 1), coord.At(i, 2))
		S.h.Write([]byte(str))
	}
	S.frames++
	return nil
}

//Only the first map will be read!
func NewWriter(name string, mol chem.Atomer, header map[string]string, compressionLevel ...int) (*CsfW, error) {
	natoms := mol.Len()
	var level int = 5
	if len(compressionLevel) > 0 {
		level = compressionLevel[0]
	}
	S := new(CsfW)
	var err error
	S.f, err = os.Create(name)
	if err != nil {
		return nil, err
	}

	//We need the symbols for the XYZ
	S.symbols = make([]string, mol.Len())
	for i := 0; i < mol.Len(); i++ {
		sym := mol.Atom(i).Symbol
		if sym == "" {
			sym = "X" //It doesn't actually matter, it's just for XYZ compatibility
		}
		S.symbols[i] = sym
	}

	format := strings.ToLower(name)[len(name)-1]
	zwriter := func(a io.Writer) (io.WriteCloser, error) {
		r, err := flate.NewWriter(a, level)
		return r, err
	}

	var AnyNewWriter func(io.Writer) (io.WriteCloser, error)
	switch format {
	case 'l':
		AnyNewWriter = func(a io.Writer) (io.WriteCloser, error) { return lzw.NewWriter(a, lzw.MSB, lzwLitwidth), nil }
	case 'f':
		AnyNewWriter = zwriter
	case 'z':
		AnyNewWriter = func(a io.Writer) (io.WriteCloser, error) { return gzip.NewWriterLevel(a, level) }
	default:
		AnyNewWriter = zwriter

	}

	S.h, err = AnyNewWriter(S.f)
	if err != nil {
		return nil, Error{"Can't read header " + err.Error(), S.filename, []string{"NewWriter"}, true}

	}

	//	S.h = AnyNewWriter(S.f)
	S.natoms = natoms
	S.filename = name
	S.writeable = true
	S.header = header

	/*
		if header != nil {
			headerstr := ""
			for k, v := range header {
				headerstr += fmt.Sprintf("%s=%v\n", k, v)
			}
			S.h.Write([]byte(headerstr))

		}
		S.h.Write([]byte(fmt.Sprintf("** %d\n", S.natoms)))
	*/
	return S, nil
}

//Read!
type CsfR struct {
	f            *os.File
	lzw          io.ReadCloser
	h            *bufio.Reader
	intermediate *bufio.Reader
	natoms       int
	filename     string
	readable     bool
	framesread   int
}

func New(name string) (*CsfR, map[string]string, error) {
	S := new(CsfR)
	S.natoms = -1 //just so we know if things don't work
	var err error
	S.filename = name
	S.f, err = os.Open(S.filename)
	if err != nil {
		return nil, nil, err
	}
	var AnyNewReader func(io.Reader) (io.ReadCloser, error)
	zreader := func(a io.Reader) (io.ReadCloser, error) {
		r := flate.NewReader(a)
		return r, nil
	}

	switch strings.ToLower(name)[len(name)-1] {
	case 'l':
		AnyNewReader = func(a io.Reader) (io.ReadCloser, error) { return lzw.NewReader(a, lzw.MSB, 8), nil }
	case 'f':
		AnyNewReader = zreader
	case 'z':
		AnyNewReader = func(a io.Reader) (io.ReadCloser, error) { return gzip.NewReader(a) }
	default:
		AnyNewReader = zreader
	}

	S.intermediate = bufio.NewReader(S.f)
	S.lzw, err = AnyNewReader(S.intermediate)
	if err != nil {
		return nil, nil, Error{"Can't read header " + err.Error(), S.filename, []string{"NewCsfR"}, true}
	}
	S.h = bufio.NewReader(S.lzw)
	str, err := S.h.ReadString('\n')
	if err != nil {
		return nil, nil, Error{"Can't read Atom number " + err.Error(), S.filename, []string{"NewCsfR"}, true}
	}
	str = strings.TrimSuffix(str, "\n")
	S.natoms, err = strconv.Atoi(str)
	if err != nil {
		return nil, nil, Error{fmt.Sprintf("Can't parse atom number from '%s': %s", str, err.Error()), S.filename, []string{"NewCsfR"}, true}
	}
	header, err := S.h.ReadString('\n')
	if err != nil {
		return nil, nil, Error{"Can't read header " + err.Error(), S.filename, []string{"NewCsfR"}, true}
	}
	header = strings.TrimSuffix(header, "\n")
	sheader := strings.Split(header, "**")
	rmap := make(map[string]string)
	if len(sheader) == 0 {

		return S, rmap, nil
	}
	for _, v := range sheader {
		s := strings.Split(v, "=")
		if len(s) > 1 {
			rmap[s[0]] = s[1]
		}
	}
	return S, rmap, nil
}

func (S *CsfR) frameError(err error) error {
	if err != nil {
		if err.Error() == "EOF" {
			//nothing bad happened here, the trajectory just ended.
			S.Close()
			return newlastFrameError(S.filename, "Next")
		}
		return Error{"Can't read frame " + err.Error(), S.filename, []string{"Next"}, true}
	}
	return nil
}

func (S *CsfR) Next(c *v3.Matrix) error {
	//the header for the first frame is read upon initialization, so we skip it here.
	if S.framesread > 0 {
		_, err := S.h.ReadBytes('\n') //the natoms line
		if err != nil {
			return S.frameError(err)
		}
		_, err = S.h.ReadBytes('\n') //the empty line
		if err != nil {
			//we shouldn't find an EOF here, so I'll return a proper error
			return Error{"Can't read frame " + err.Error(), S.filename, []string{"Next"}, true}
		}

	}
	for i := 0; i < S.natoms; i++ {
		b, err := S.h.ReadBytes('\n')
		if err != nil {
			if err.Error() == "EOF" {
				//nothing bad happened here, the trajectory just ended.
				S.Close()
				return newlastFrameError(S.filename, "Next")
			}
			return Error{"Can't read frame " + err.Error(), S.filename, []string{"Next"}, true}
		}
		str := (b[0 : len(b)-1]) //This should work for ASCII, which is fine for Csf. It removes the last '\n'
		coords := bytes.Fields(str)
		if len(coords) != 4 { //This also checks for a premature "end of frame" line, i.e. "*", which would have len 1
			return Error{"Wrong number of coordinates or atoms in frame: " + string(str), S.filename, []string{"Next"}, true}
		}
		if c == nil {
			continue //We ignore this whole frame, reading the content but not saving it.
			//Note that we still check the frame for correctness.
		}
		for j, v := range coords[1:] {
			//This should be replaced for a call to bytesconv.BytesToString when I update the compiler.
			n, err := strconv.ParseFloat(*(*string)(unsafe.Pointer(&v)), 64)
			if err != nil {
				return Error{"Coordinate un parseable in frame." + err.Error(), S.filename, []string{"Next"}, true}
			}
			c.Set(i, j, n)
		}

	}
	S.framesread++
	return nil
}

func (S *CsfR) Close() {
	if !S.readable {
		return
	}
	S.lzw.Close()
	S.readable = false
	return
}

func (S *CsfR) Len() int {
	return S.natoms
}

//Errors

//Errors

//errDecorate is a helper function that asserts that the error is
//implements chem.Error and decorates the error with the caller's name before returning it.
//if used with a non-chem.Error error, it will cause a panic.
func errDecorate(err error, caller string) error {
	err2 := err.(chem.Error) //I know that is the type returned byt initRead
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
	return fmt.Sprintf("stf file %s error: %s", err.filename, err.message)
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
	TrajUnIniRead       = "Traj object uninitialized to read"
	TrajUnIniWrite      = "Traj object uninitialized to read"
	ReadError           = "Error reading frame"
	UnableToOpen        = "Unable to open file"
	SecurityCheckFailed = "Failed Security Check"
	NilCoordinates      = "Given nil coordinates"
	WrongFormat         = "Wrong format in the STF file or frame"
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
