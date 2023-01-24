package histo

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"sort"
	"strings"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/stat"
)

//Combines 2 matrices element-wise using the function f, which should take 2 histograms  to be
//combined and one more where the result of the operation is stored.
func MatrixCombine(f func(a, b, dest *Data), a, b, dest *Matrix) {
	if a.rows != b.rows || a.cols != b.cols || a.rows != dest.rows || a.cols != dest.cols {
		panic("goChem/histo.MatrixMerge: Ill-formed matrices for merging")
	}
	//This should work if they are both nil
	if !(a.dividers == nil && b.dividers == nil) && !floats.Equal(a.dividers, b.dividers) {
		panic("goChem/histo.MatrixMerge: Matrices don't have the same dividers")
	}
	for i, v := range dest.d {
		f(a.d[i], b.d[i], v)
	}
}

//A matrix of histograms
type Matrix struct {
	rows, cols int       //total
	d          []*Data   //row-major
	dividers   []float64 //if not nil, all histograms have the same dividers
}

//NewMatrix returns a new matrix of *Data with r and c rows and column
//and dividers dividers. Dividers can be nil, in which case, elements
//of the matrix will not be forced to have the same dividers
func NewMatrix(r, c int, dividers []float64) *Matrix {
	ret := new(Matrix)
	ret.rows = r
	ret.cols = c
	ret.d = make([]*Data, r*c)
	ret.dividers = dividers
	return ret
}

func (M *Matrix) Dims() (int, int) {
	return M.rows, M.cols
}

//Copies the dividers of the histogram
func (M *Matrix) CopyDividers(dest ...[]float64) []float64 {
	if M.dividers == nil {
		return nil
	}
	d := getCopySlice(len(M.dividers), dest...)
	return floats.ScaleTo(d, 0, M.dividers)
}

func (M *Matrix) String() string {
	ret := fmt.Sprintf("rows:%d cols:%d | Data:\n", M.rows, M.cols)
	t := make([]string, 0, len(M.d))
	for _, v := range M.d {
		t = append(t, v.String())
	}
	return ret + strings.Join(t, "\n\n")
}

func (M *Matrix) MarshalJSON() ([]byte, error) {
	j, err := json.Marshal(struct {
		Rows     int       `json:"rows"`
		Cols     int       `json:"cols"`
		D        []*Data   `json:"data"`
		Dividers []float64 `json:"dividers"`
	}{
		Rows:     M.rows,
		Cols:     M.cols,
		D:        M.d,
		Dividers: M.dividers,
	})
	if err != nil {
		return nil, err
	}
	return j, nil
}

func (M *Matrix) UnmarshalJSON(b []byte) error {
	var a struct {
		Rows     int       `json:"rows"`
		Cols     int       `json:"cols"`
		D        []*Data   `json:"data"`
		Dividers []float64 `json:"dividers"`
	}

	err := json.Unmarshal(b, &a)
	if err != nil {
		return err
	}
	M.rows = a.Rows
	M.cols = a.Cols
	M.d = a.D
	M.dividers = a.Dividers
	return nil
}

//returns the index in the []*Data slice of a matrix given
//the row and column indexes.
//just to avoid fixing it in many places if I screw up
func (M *Matrix) rc2i(r, c int) int {
	M.Check(r, c, true)
	return M.cols*r + c
}

//Fill fills the matrix with empty histograms
//If the matrix has a non-nil delimiters slice,
//that slice is used for all the histograms created
func (M *Matrix) Fill() {
	for i := 0; i < M.rows; i++ {
		for j := 0; j < M.cols; j++ {
			M.NewHisto(i, j, M.dividers, nil)
		}

	}
}

//Check checks if the given row and column indexes are within range.
//if pan is given and true, it panics if either is out of range,
//otherwise, it returns an error.
func (M *Matrix) Check(r, c int, pan ...bool) error {
	p := false
	var err error
	if len(pan) > 0 && pan[0] {
		p = true
	}
	if r >= M.rows {
		err = fmt.Errorf("goChem/Histo: Row out of range")
	}
	if c >= M.cols {
		err = fmt.Errorf("goChem/Histo: Column out of range")
	}
	if err != nil && p {
		panic(err.Error())
	}
	return err
}

//NewHisto Puts a new histogram in the r,c position in the matrix. Dividers can be nil, in which case, the matrix
//should have its dividers. If there are no dividers, or they don't match, the function will panic.
//rawdata can also be nil, in which case, an empty histogram will be put in the position.
func (M *Matrix) NewHisto(r, c int, dividers []float64, rawdata []float64, ID ...int) {
	if dividers == nil {
		if M.dividers != nil {
			dividers = M.dividers
		} else {
			panic("goChem/histo.Matrix.NewHisto: dividers not given, and can't be taken from current first element")
		}
	} else if M.dividers != nil && !floats.Equal(M.dividers, dividers) {
		//Maybe this should be returned as an error instead
		log.Printf("goChem/histo.Matrix.NewHisto: dividers given but don't match the dividers of the matrix. The matrix's dividers will be used.")
		dividers = M.dividers
	}
	M.d[M.rc2i(r, c)] = NewData(dividers, rawdata,ID...)
}

//View Returns a view of the histogram in the r,c position in the matrix
func (M *Matrix) View(r, c int) *Data {
	return M.d[M.rc2i(r, c)]
}

//Adds one or more data points to the histogram in the r,c position in the matrix
func (M *Matrix) AddData(r, c int, point ...float64) {
	M.d[M.rc2i(r, c)].AddData(point...)
}

//Normalize all the histograms in the matrix
func (M *Matrix) NormalizeAll() {
	for _, v := range M.d {
		v.Normalize()
	}
}

//Un-normalize all the histograms in the matrix
func (M *Matrix) UnNormalizeAll() {
	for _, v := range M.d {
		v.UnNormalize()
	}
}

//Applies the f function to each element in the matrix, the results are returned as
//a [][]float64. Also returns error unpon failure, or nil.
func (M *Matrix) FromAll(f func(D *Data) (float64, error)) ([][]float64, error) {
	r := make([][]float64, M.rows)
	var err error
	for i := 0; i < M.rows; i++ {
		r[i] = make([]float64, M.cols)
		for j := 0; j < M.cols; j++ {
			r[i][j], err = f(M.d[M.rc2i(i, j)])
			if err != nil {
				return nil, fmt.Errorf("goChem/Histo.Matrix.FromAll: Error at %d, %d: %v", i, j, err)
			}
		}
	}
	return r, nil
}

//Applies the f function to each element in the matrix. Returns error unpon failure, or nil.
func (M *Matrix) ToAll(f func(D *Data) error) error {
	var err error
	for i := 0; i < M.rows; i++ {
		for j := 0; j < M.cols; j++ {
			err = f(M.d[M.rc2i(i, j)])
			if err != nil {
				return fmt.Errorf("goChem/Histo.Matrix.ToAll: Error at %d, %d: %v", i, j, err)
			}

		}
	}
	return nil
}

type Data struct {
    	id int
	normalized bool
	total      int
	dividers   []float64
	histo      []float64
}

func (D *Data) MarshalJSON() ([]byte, error) {
	j, err := json.Marshal(struct {
	    	ID         int       `json:"id"`
		Normalized bool      `json:"normalized"`
		Total      int       `json:"total"`
		Dividers   []float64 `json:"dividers"`
		Histo      []float64 `json:"histo"`
	}{
	    	ID: D.id
		Normalized: D.normalized,
		Total:      D.total,
		Dividers:   D.dividers,
		Histo:      D.histo,
	})
	if err != nil {
		return nil, err
	}
	return j, nil
}

func (D *Data) UnmarshalJSON(b []byte) error {
	var a struct {
	    	ID         int       `json:"id"`
		Normalized bool      `json:"normalized"`
		Total      int       `json:"total"`
		Dividers   []float64 `json:"dividers"`
		Histo      []float64 `json:"histo"`
	}

	err := json.Unmarshal(b, &a)
	if err != nil {
		return err
	}
	D.id         = a.ID
	D.normalized = a.Normalized
	D.total = a.Total
	D.dividers = a.Dividers
	D.histo = a.Histo
	return nil
}


//ID returns the ID of the histogram
func (D *Data) ID() int {
	return D.id
}


//String prints a -hopefully- pretty string representation of
//the histogram. The representation uses 3 lines of thext
func (D *Data) String() string {
    ret := fmt.Sprintf("ID: %d, Normalized: %v, TotalData: %d\n", D.id, D.normalized, D.total)
	d := make([]string, 0, len(D.dividers)-1)
	h := make([]string, 0, len(D.dividers)-1)
	for i, v := range D.histo {
		d = append(d, fmt.Sprintf("%4.2f-%4.2f", D.dividers[i], D.dividers[i+1]))
		h = append(h, fmt.Sprintf("%9.3f", v))
	}
	//fmt.Println(h, D.histo) /////////
	return ret + fmt.Sprintf("%s\n%s", strings.Join(d, " "), strings.Join(h, " "))

}

//Returns a new histogram from the dividers and rawdata given
//rawdata can be nil. In that case, an empty histogram is created.
//if an ID for the histogram is given, it will be set. If not, the ID will
//be set to -1.
func NewData(dividers []float64, rawdata []float64, ID ...int) *Data {
	d := new(Data)
	//I prefer to copy the slice to avoid somebody changing it from outside
	d.dividers = make([]float64, len(dividers))
	for i, v := range dividers {
		d.dividers[i] = v
	}
	d.histo = make([]float64, len(dividers)-1)
	if rawdata != nil {
		//d.total = len(rawdata)
		d.ReHisto(d.dividers, rawdata)
		//println("not nil!") ///////
	}
	d.id=-1
	if len(ID)>0{
		d.id=ID[0]
	}
	return d

}

//Adds the given data point(s) to the histogram
func (M *Data) AddData(point ...float64) {
	var norma bool
	if M.normalized {
		norma = true
		M.UnNormalize()
	}
	for _, v := range point {
		for j, w := range M.dividers {
			//Values that are larger than the last divider are just omitted.
			if j == len(M.dividers)-1 {
				break
			}
			if w <= v && v < M.dividers[j+1] {
				M.histo[j]++
				break
			}
		}
	}
	M.total += len(point)
	//if it was normalized, we should return it to that state
	if norma {
		M.Normalize()
	}
}

//Normalized Returns true if the histogram is normalized
func (D *Data) Normalized() bool {
	return D.normalized
}

//Normalize normalizes the histogram
func (D *Data) Normalize() {
	D.normaunnorma(true)
}

//UnNormalize un-normalizes the histogram
func (D *Data) UnNormalize() {
	D.normaunnorma(false)
}

//normalizes or un-normalizes the histogram depending
//on whether normalize is true
func (D *Data) normaunnorma(normalize bool) {
	if D.total <= 0 {
		return
	}
	n := float64(D.total)
	D.normalized = false
	if normalize {
		n = 1 / float64(D.total)
		D.normalized = true
	}

	floats.Scale(n, D.histo)

}

//Copies the dividers of the histogram
func (D *Data) CopyDividers(dest ...[]float64) []float64 {
	d := getCopySlice(len(D.dividers), dest...)
	return floats.ScaleTo(d, 0, D.dividers)
}

func (D *Data) Copy(dest ...[]float64) []float64 {
	d := getCopySlice(len(D.histo), dest...)
	return floats.ScaleTo(d, 0, D.histo)
}

func (D *Data) View() []float64 {
	return D.histo
}

//Add adds the histograms a and b putting the result in the receiver.
func (D *Data) Add(a, b *Data) {
	D.dividers = a.CopyDividers(D.dividers)
	if len(a.dividers) != len(b.dividers) {
		panic("goChem/Histo.Data.Add: Ill-formed histograms for addition")
	}

	for i, v := range a.dividers {
		if v != b.dividers[i] {
			panic("goChem/Histo.Data.Add: Dividers must match in added histograms")
		}
		if i == len(a.dividers)-1 {
			break //a.histo has 1 less element than a.dividers, so we skip the next operation for the last one.
		}
		D.histo[i] = a.histo[i] + b.histo[i]
	}
}

//Sub substract the histograms a and b puting the results in the receiver
//if abs is given and true (only the first element is considered)
func (D *Data) Sub(a, b *Data, abs ...bool) {
	f := func(a float64) float64 { return a }
	if len(abs) > 0 && abs[0] {
		f = func(a float64) float64 { return math.Abs(a) }
	}
	D.dividers = a.CopyDividers(D.dividers)
	if len(a.dividers) != len(b.dividers) {
		panic("goChem/Histo.Data.Add: Ill-formed histograms for addition")
	}

	for i, v := range a.dividers {
		if v != b.dividers[i] {
			panic("goChem/Histo.Data.Add: Dividers must match in added histograms")
		}
		if i == len(a.dividers)-1 {
			break //a.histo has 1 less element than a.dividers, so we skip the next operation for the last one.
		}

		D.histo[i] = f(a.histo[i] - b.histo[i])
	}
}

func (D *Data) Sum() float64 {
	return floats.Sum(D.histo)
}

func (D *Data) ReHisto(dividers, rawdata []float64) {
	if rawdata != nil {
		sort.Float64s(rawdata)
		//stat.Histograms just panics instead of omitting the values that are off limits
		//so we remove them here before the call.
		maxi := sort.SearchFloat64s(rawdata, dividers[len(dividers)-1])
		mini := sort.SearchFloat64s(rawdata, dividers[0])
		if maxi < len(rawdata) {
			rawdata = rawdata[:maxi]
		}
		if mini != 0 {
			rawdata = rawdata[mini:]
		}

	}
	D.total = len(rawdata) //as this could have been modified
	D.histo = stat.Histogram(nil, dividers, rawdata, nil)
	//println(D.histo[0]) ////////////
}

func getCopySlice(N int, dest ...[]float64) []float64 {
	var d []float64
	if len(dest) > 0 && len(dest[0]) >= N {
		d = dest[0]
		if len(dest[0]) > N {
			d = dest[0][:N] //floats.ScaleTo wants both slices to _match_
		}
	} else {
		d = make([]float64, N)
	}
	return d

}
