package stf

import (
	"bufio"
	"compress/flate"
	"compress/gzip"
	"compress/lzw"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
	"unsafe"

	"github.com/klauspost/compress/zstd"
	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/mat"
)

const (
	lzwLitwidth int = 8
)

//Write!
type StfW struct {
	f           *os.File
	h           io.WriteCloser
	natoms      int
	filename    string
	writeable   bool
	framebuffer *v3.Matrix
	prec        int
}

func (S *StfW) Close() {

	if S == nil {
		return
	}
	if S.writeable {
		S.h.Close()
		S.f.Close()
	}
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

func (S *StfW) WNext(coord *v3.Matrix, box ...[]float64) error {
	//centroid, err := chem.MassCenterMem(coord, coord, S.framebuffer) //not actually the CoM, but the geometric center.
	//if err == nil {
	//	coord = S.framebuffer //we won't say anything in case of error, sorry.
	//}

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
	//	strs := make([]string, 4) //old code
	var temp [3]int
	var floats [3]float64
	var str string
	//	prec = 2 //default
	for i := 0; i < v; i++ {
		floats[0] = coord.At(i, 0)
		floats[1] = coord.At(i, 1)
		floats[2] = coord.At(i, 2)

		str = coordsEncode(floats, temp, S.prec)

		S.h.Write([]byte(str))
	}
	if len(box) > 0 && len(box[0]) >= 9 {
		b := box[0]
		//if we did do the centroid thing, we should also displace the box vectors.
	//	if centroid != nil {
	//		for in := 0; in < 3; in++ {
	//			b[(3*in)+0] -= centroid.At(0, 0) //I like it better with the parenthesis :-D
	//			b[(3*in)+1] -= centroid.At(0, 1)
	//			b[(3*in)+2] -= centroid.At(0, 2)
//
//			}
		}
		S.h.Write([]byte(fmt.Sprintf("* %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n", b[0],
			b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8])))

	} else {
		S.h.Write([]byte("*\n"))
	}
	//	println("Wrote a frame!") ///////////////////////////
	return nil
}

//Only the first map will be read!
func NewWriter(name string, natoms int, header map[string]string, compressionLevel ...int) (*StfW, error) {
	var level int = 11 //For python compatibility
	if len(compressionLevel) > 0 {
		level = compressionLevel[0]
	}
	S := new(StfW)
	var err error
	S.f, err = os.Create(name)
	if err != nil {
		return nil, err
	}
	S.framebuffer = v3.Zeros(natoms)
	format := strings.ToLower(name)[len(name)-1]
	zwriter := func(a io.Writer) (io.WriteCloser, error) {
		r, err := flate.NewWriter(a, level)
		return r, err
	}
	gzipwriter := func(a io.Writer) (io.WriteCloser, error) { return gzip.NewWriterLevel(a, level) }
	zstdwriter := func(a io.Writer) (io.WriteCloser, error) {
		return zstd.NewWriter(a, zstd.WithEncoderLevel(zstd.SpeedBestCompression))
	}
	//	snappywriter := func(a io.Writer) (io.WriteCloser, error) { return snappy.NewBufferedWriter(a), nil }

	var AnyNewWriter func(io.Writer) (io.WriteCloser, error)
	switch format {
	case 'l':
		AnyNewWriter = func(a io.Writer) (io.WriteCloser, error) { return lzw.NewWriter(a, lzw.MSB, lzwLitwidth), nil }
	case 'f':
		AnyNewWriter = zstdwriter
	case 'z':
		AnyNewWriter = gzipwriter
	case 's':
		AnyNewWriter = zstdwriter
	case 'r':
		AnyNewWriter = zwriter

	default:
		AnyNewWriter = zstdwriter

	}

	S.h, err = AnyNewWriter(S.f)
	if err != nil {
		return nil, Error{"Can't read header " + err.Error(), S.filename, []string{"NewWriter"}, true}

	}

	//	S.h = AnyNewWriter(S.f)
	S.natoms = natoms
	S.filename = name
	S.writeable = true
	S.prec = 2 //the default
	if header != nil {
		if p, ok := header["prec"]; ok && p != "2" {
			prec, err := strconv.Atoi(p)
			if err != nil {
				S.prec = prec
			} else {
				log.Printf("Invalid precision for trajectory %s. Will use the default", S.filename)
			}
		}
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
	//	framebuffer  *v3.Matrix
	natoms   int
	filename string
	prec     int
	readable bool
}

//This will cause additional indirections
//but I suppose it won't matter, as each call will
//take enough time to make those delays irrelevant.
//Also, why couldn't *zstd.Decoder implement io.ReadCloser? :-(
type stdql struct {
	closeql func() //The things I have to do xD
	*zstd.Decoder
}

//Close Closes the object. It can not be used after this call
func (s stdql) Close() error {
	s.closeql()
	return nil
}

func coordsEncode(f [3]float64, temp [3]int, prec int) string {
	p := 100.0
	if prec > 0 && prec != 2 { //2 is the current value, so we do nothign in that case
		p = math.Pow(10.0, float64(prec))
	}
	for i, v := range f {
		temp[i] = int(math.RoundToEven(v * p))
	}
	return fmt.Sprintf("%d %d %d\n", temp[0], temp[1], temp[2])
}

/*
func hexaDecEncode(f [3]float64, temp [3]int, ts [3][2]string, prec int, big ...bool) string {
	p := 100.0
	var ret string
	if prec > 0 && prec != 2 { //2 is the current value, so we do nothign in that case
		p = math.Pow(10.0, float64(prec))
	}
	for i, v := range f {
		temp[i] = int(math.RoundToEven(v * p))
	}
	if len(big) > 0 && big[0] {
		ret = fmt.Sprintf("%x %x %x", temp[0], temp[1], temp[2])
		//All this mess is to ensure that the numbers use either 4 spaces, if positive, or 5, if negative, but never less.
	} else {
		for i, v := range temp {
			if v < 0 {
				ts[i][0] = fmt.Sprintf("%04x", v*-1)
				ts[i][1] = "-"
			} else {
				ts[i][0] = fmt.Sprintf("%04x", v)
				ts[i][1] = ""
			}
		}
		ret = fmt.Sprintf("%s%s%s%s%s%s\n", ts[0][1], ts[0][0], ts[1][1], ts[1][0], ts[2][1], ts[2][0])
	}
	//	for i, _ := range temp {
	//		temp[i] = 0 //I really don't want to risk any contamination here
	//	}
	return ret

}




func mustParseHexInt(s string) int64 {
	ret, err := strconv.ParseInt(strings.Replace(s, " ", "", -1), 16, 32)
	if err != nil {
		panic(err.Error()) //really shouldn't happen
	}
	return ret
}


func hexaDecDecode(str string, temp [3]float64, coords [3]string, delays [3]int, prec int, big ...bool) {
	mult := 100.0
	if prec > 0 && prec != 2 { //2 is just the current value, so we can save the operation
		mult = math.Pow(10, float64(prec))
	}
	delays[0] = 0 //just to be sure
	delays[1] = 0
	delays[2] = 0
	str = strings.Replace(str, "\n", "", -1) //maybe wasteful, but I'd rather be careful. We can see if things get too slow
	if len(big) > 0 && big[0] {
		s := strings.Fields(str)
		for i, v := range s {
			coords[i] = v
		}
		//the most common case. This
		//whole thing is to deal with the fact that a "-" sign means
		//that the coresponding field will take 5 spaces instead of 4.
	} else {
		if str[0] == '-' {
			delays[0] = 1
		} else {

		}
		coords[0] = str[0 : 5+delays[0]]
		if str[5+delays[0]] == '-' {
			delays[1] = 1
		}
		coords[1] = str[5+delays[0] : 9+delays[0]+delays[1]]
		if str[9+delays[0]+delays[1]] == '-' {
			delays[2] = 1
		}
		coords[2] = str[9+delays[0]+delays[1]:]
	}
	for i, v := range coords {
		temp[i] = float64(mustParseHexInt(v)) / mult

	}

}
**************/
//New opens a STF trajectory for reading, and returns a pointer
//to the handle, a map with the metadata (or nil, if no metadata is found)
//and error or nil.
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
	zstdreader := func(a io.Reader) (io.ReadCloser, error) {
		r, err := zstd.NewReader(a)
		var ql *stdql
		ql = &stdql{r.Close, r}

		return ql, err

	}
	gzreader := func(a io.Reader) (io.ReadCloser, error) { return gzip.NewReader(a) }
	switch strings.ToLower(name)[len(name)-1] {
	case 'l':
		AnyNewReader = func(a io.Reader) (io.ReadCloser, error) { return lzw.NewReader(a, lzw.MSB, 8), nil }
	case 'f':
		AnyNewReader = zstdreader
	case 'z':
		AnyNewReader = gzreader
	case 's':
		AnyNewReader = zstdreader
	case 'r':
		AnyNewReader = zreader

	default:
		AnyNewReader = zstdreader

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
	//	S.framebuffer = v3.Zeros(S.natoms)
	S.readable = true
	if m != nil {
		if p, ok := m["prec"]; ok && p != "2" {
			prec, err := strconv.Atoi(p)
			if err != nil {
				S.prec = prec
			} else {
				log.Printf("Invalid precision for trajectory %s. Will assume the default", S.filename)
			}
		}

	}
	return S, m, nil
}

//Readabe returns true if the handle is readable (if it is possible to call Next on it)
func (S *StfR) Readable() bool {
	return S.readable
}

func coordsDecode(str string, temp *[3]float64, prec int) error {
	p := 100.0
	if prec > 0 && prec != 2 { //2 is just the current value, so we can save the operation
		p = math.Pow(10.0, float64(prec))

	}
	s := strings.Fields(str)
	if len(s) < 3 {
		return fmt.Errorf("Ill formated coordinates line in stf: Too few fields: %s", str)
	}
	if len(s) > 3 {
		return fmt.Errorf("Ill formated coordinates line in stf: Too many fields: %s", str)
	}
	for i, v := range s {
		f, err := strconv.Atoi(v)
		if err != nil {
			return fmt.Errorf("Can't parse coordinate %d (%s). Error: %s", i, v, err.Error())
		}
		temp[i] = float64(f) / p
	}
	return nil
}

//Next puts in the given matrix (c) the coordinates for the next frame of the trajectory
//and, if given, and the information is present, puts the box vector information in box
//Returns error if the operation is not successful. If the error is EOF, the end of the
//trajectory has been reached, not an actual error.
func (S *StfR) Next(c *v3.Matrix, box ...[]float64) error {
	//	var tmpstr string
	//	var prec int = 2
	//	var pointplace int
	var temp [3]float64
	var err error
	for i := 0; i < S.natoms; i++ {
		b, err := S.h.ReadBytes('\n')
		if err != nil {
			// EOF should only happen when reading the first atom
			if err.Error() == "EOF" && (strings.Contains(err.Error(), "EOF")) && i == 0 {
				//nothing bad happened here, the trajectory just ended.
				S.Close()
				return newlastFrameError(S.filename, "Next")
			} else {
				return Error{message: err.Error(), filename: S.filename, critical: true}
			}

		}
		err = coordsDecode(string(b[:len(b)-1]), &temp, S.prec)
		if err != nil {
			return Error{message: err.Error(), filename: S.filename, critical: true}
		}

		if c == nil {
			continue //We ignore this whole frame, reading the content but not saving it.
			//Note that we still check the frame for correctness.
		}
		for j, v := range temp {
			c.Set(i, j, v)
		}

	}
	s, err := S.h.ReadString('\n')
	if err != nil {
		return Error{"Can't read the frame termination mark" + err.Error(), S.filename, []string{"Next"}, true}
	}
	if s[0] != '*' {
		return Error{"Wrong number of atoms in frame" + err.Error(), S.filename, []string{"Next"}, true}
	}

	//This part reads the
	if len(box) > 0 && len(box[0]) >= 9 {
		fields := strings.Fields(strings.TrimSpace(s)) //we get rid of that last "\n"
		if len(fields) >= 10 {                         // The "*" and the 9 numbers
			var errbox error
			for j, v := range fields[1:] {
				//This should be replaced for a call to bytesconv.BytesToString when I update the compiler.
				box[0][j], errbox = strconv.ParseFloat(*(*string)(unsafe.Pointer(&v)), 64)
				if errbox != nil {
					break
				}
			}
			//If we got an error reading any of the values, we just set the whole thing to zero
			//and log, no error returned.
			if errbox != nil {
				log.Printf("Failed to read box in a frame from %s", S.filename) //just a head-up
				for i, _ := range box[0] {
					box[0][i] = 0.0
				}
			}

		} else {
			log.Printf("Trajectory file %s does not contain (correct) box information: %s", S.filename, fields) //just a head-up
		}
	}

	return nil
}

//Close closes the object, and marks it as unreadable
func (S *StfR) Close() {
	if !S.readable {
		return
	}
	S.lzw.Close()
	S.readable = false
	return
}

//Len returns the number of atoms in each frame of the trajectory.
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
