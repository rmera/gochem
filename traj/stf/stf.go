package stf

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

const (
	lzwLitwidth int = 8
)

//Write!
type StfW struct {
	f         *os.File
	h         io.WriteCloser
	natoms    int
	filename  string
	writeable bool
}

func (S *StfW) Close() {
	S.h.Close()
	S.f.Close()
	S.writeable = false
	return
}

func (S *StfW) Len() int {
	return S.natoms
}

//compatibility with Gonum
func (S *StfW) WNextDense(dcoord *mat.Dense) error {
	coord := v3.Dense2Matrix(dcoord)
	err := S.WNext(coord)
	if err != nil {
		err = errDecorate(err, "WNextDense")
	}
	return err
}

func (S *StfW) WNext(coord *v3.Matrix) error {
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
	strs := make([]string, 4)
	for i := 0; i < v; i++ {
		ff := strconv.FormatFloat
		prec := 3
		strs[0] = ff(coord.At(i, 0), 'f', prec, 64)
		strs[1] = ff(coord.At(i, 1), 'f', prec, 64)
		strs[2] = ff(coord.At(i, 2), 'f', prec, 64)
		strs[3] = "\n"
		str := strings.Join(strs, " ")
		//str := fmt.Sprintf("%07.3f %07.3f %07.3f\n", coord.At(i, 0), coord.At(i, 1), coord.At(i, 2))
		S.h.Write([]byte(str))
	}
	S.h.Write([]byte("*\n"))
	return nil
}

//Only the first map will be read!
func NewWriter(name string, natoms int, header map[string]string, compressionLevel ...int) (*StfW, error) {
	var level int = 5 //For python compatibility :-/
	if len(compressionLevel) > 0 {
		level = compressionLevel[0]
	}
	S := new(StfW)
	var err error
	S.f, err = os.Create(name)
	if err != nil {
		return nil, err
	}

	format := strings.ToLower(name)[len(name)-1]
	zwriter := func(a io.Writer) (io.WriteCloser, error) {
		r, err := flate.NewWriter(a, level)
		return r, err
	}
	gzipwriter := func(a io.Writer) (io.WriteCloser, error) { return gzip.NewWriterLevel(a, level) }

	var AnyNewWriter func(io.Writer) (io.WriteCloser, error)
	switch format {
	case 'l':
		AnyNewWriter = func(a io.Writer) (io.WriteCloser, error) { return lzw.NewWriter(a, lzw.MSB, lzwLitwidth), nil }
	case 'f':
		AnyNewWriter = zwriter
	case 'z':
		AnyNewWriter = gzipwriter
	default:
		AnyNewWriter = gzipwriter

	}

	S.h, err = AnyNewWriter(S.f)
	if err != nil {
		return nil, Error{"Can't read header " + err.Error(), S.filename, []string{"NewWriter"}, true}

	}

	//	S.h = AnyNewWriter(S.f)
	S.natoms = natoms
	S.filename = name
	S.writeable = true
	if header != nil {
		headerstr := ""
		for k, v := range header {
			headerstr += fmt.Sprintf("%s=%v\n", k, v)
		}
		S.h.Write([]byte(headerstr))

	}
	S.h.Write([]byte(fmt.Sprintf("** %d\n", S.natoms)))
	return S, nil
}

//Read!
type StfR struct {
	f            *os.File
	lzw          io.ReadCloser
	h            *bufio.Reader
	intermediate *bufio.Reader
	natoms       int
	filename     string
	readable     bool
}

func New(name string) (*StfR, map[string]string, error) {
	S := new(StfR)
	S.natoms = -1 //just so we know if things don't work
	var m map[string]string
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
		return nil, nil, Error{"Can't read header " + err.Error(), S.filename, []string{"NewStfR"}, true}
	}
	S.h = bufio.NewReader(S.lzw)
	for {
		str, err := S.h.ReadString('\n')
		if err != nil {
			return nil, nil, Error{"Can't read header " + err.Error(), S.filename, []string{"NewStfR"}, true}
		}
		str = strings.TrimSuffix(str, "\n")
		//	str = string([]byte(str)[0 : len(str)-1]) //This should work for ASCII, which is fine for Stf. It removes the last '\n'
		if strings.Contains(str, "**") {
			nat := strings.Fields(str)
			if len(nat) < 2 {
				return nil, nil, Error{fmt.Sprintf("Can't read atom number from '%s': %s", str, err.Error()), S.filename, []string{"NewStfR"}, true}
			}
			nat[1] = strings.TrimSuffix(nat[1], "\n")
			S.natoms, err = strconv.Atoi(nat[1])
			if err != nil {
				return nil, nil, Error{fmt.Sprintf("Can't read atom number from '%s': %s", nat[1], err.Error()), S.filename, []string{"NewStfR"}, true}

			}
			break
		}
		kv := strings.Split(str, "=")
		if len(kv) != 2 {
			return nil, nil, Error{"Malformed header" + err.Error(), S.filename, []string{"NewStfR"}, true}
		}
		m[kv[0]] = kv[1]
	}
	S.readable = true
	return S, m, nil
}

func (S *StfR) Readable() bool {
	return S.readable
}

func (S *StfR) Next(c *v3.Matrix) error {
	for i := 0; i < S.natoms; i++ {
		b, err := S.h.ReadBytes('\n')
		if err != nil {
			if err.Error() == "EOF" && (strings.Contains(err.Error(), "EOF") && i == 0) {
				//nothing bad happened here, the trajectory just ended.
				S.Close()
				return newlastFrameError(S.filename, "Next")
			}
			fmt.Println("i", i, b, err.Error())
			return Error{fmt.Sprintf("Can't read frame, atom: %d  read: %s: %s", i, string(b), err.Error()), S.filename, []string{"Next"}, true}
		}
		str := (b[0 : len(b)-1]) //This should work for ASCII, which is fine for Stf. It removes the last '\n'
		coords := bytes.Fields(str)
		if len(coords) != 3 { //This also checks for a premature "end of frame" line, i.e. "*", which would have len 1
			return Error{"Wrong number of coordinates or atoms in frame" + err.Error(), S.filename, []string{"Next"}, true}
		}
		if c == nil {
			continue //We ignore this whole frame, reading the content but not saving it.
			//Note that we still check the frame for correctness.
		}
		for j, v := range coords {
			//This should be replaced for a call to bytesconv.BytesToString when I update the compiler.
			n, err := strconv.ParseFloat(*(*string)(unsafe.Pointer(&v)), 64)
			if err != nil {
				return Error{"Coordinate un parseable in frame." + err.Error(), S.filename, []string{"Next"}, true}
			}
			c.Set(i, j, n)
		}

	}

	s, err := S.h.ReadString('\n')
	if err != nil {
		return Error{"Can't read the frame termination mark" + err.Error(), S.filename, []string{"Next"}, true}
	}
	if s != "*\n" {
		return Error{"Wrong number of atoms in frame" + err.Error(), S.filename, []string{"Next"}, true}
	}

	return nil
}

func (S *StfR) Close() {
	if !S.readable {
		return
	}
	S.lzw.Close()
	S.readable = false
	return
}

func (S *StfR) Len() int {
	return S.natoms
}

//NextConc takes a slice of bools and reads as many frames as elements the list has
//form the trajectory. The frames are discarted if the corresponding elemetn of the slice
//is false. The function returns a slice of channels through each of each of which
// a *matrix.DenseMatrix will be transmited
func (D *StfR) NextConc(frames []*v3.Matrix) ([]chan *v3.Matrix, error) {
	if !D.Readable() {
		return nil, Error{TrajUnIniRead, D.filename, []string{"NextConc"}, true}
	}
	framechans := make([]chan *v3.Matrix, len(frames)) //the slice of chans that will be returned
	for key, v := range frames {
		if err := D.Next(v); err != nil {
			return nil, errDecorate(err, "NextConc")
		}
		framechans[key] = make(chan *v3.Matrix)
		go func(keep *v3.Matrix, pipe chan *v3.Matrix) {
			pipe <- keep
		}(v, framechans[key])
	}
	return framechans, nil
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
