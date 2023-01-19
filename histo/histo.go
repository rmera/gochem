package histo

import (
	"fmt"
	"log"
	"math"

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
	if a.dividers != b.dividers {
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

//Returns a new matrix of *Data with r and c rows and column
//and dividers dividers. Dividers can be nil, in which case, elements
//of the matrix will not be forced to have the same dividers
func NewMatrix(r, c int, dividers []float64) *Matrix {
	ret = new(*Matrix)
	ret.rows = r
	ret.cols = c
	ret.d = make([]*Data, r*c)
	ret.dividers = dividers
	return ret
}

//returns the index in the []*Data slice of a matrix given
//the row and column indexes.
//just to avoid fixing it in many places if I screw up
func (M *Matrix) rc2i(r, c int) int {
	M.Check(r, c, true)
	return M.cols*r + c
}

//Checks if the given row and column indexes are within range.
//if pan is given and true, it panics if either is out of range,
//otherwise, it returns an error.
func (M *Matrix) Check(r, c int, pan ...bool) error {
	p := false
	var err error
	if len(pan) > 0 && pan[0] == true {
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

//Returns a view of the histogram in the r,c position in the matrix
func (M *Matrix) NewHisto(r, c int, dividers []float64, rawdata [][]float64) {
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
	M.d[M.rc2i(r, c)] = NewData(dividers, rawdata)
}

//Returns a view of the histogram in the r,c position in the matrix
func (M *Matrix) View(r, c int) *Data {
	return M.d[M.rc2i(r, c)].View()
}

//Adds one or more data points to the histogram in the r,c position in the matrix
func (M *Matrix) AddData(r, c int, point ...float64) *Data {
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
	var r [][]float64
	r = make([][]float64, M.rows)
	var err error
	for i := 0; i < M.rows; i++ {
		r[i] = appen(r[i], make([]float64, cols))
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
				return nil, fmt.Errorf("goChem/Histo.Matrix.ToAll: Error at %d, %d: %v", i, j, err)
			}

		}
	}
}

type Data struct {
	normalized bool
	total      int
	dividers   []float64
	histo      []float64
}

//Returns a new histogram from the dividers and rawdata given
//rawdata can be nil. In that case, an empty histogram is created.
func NewData(dividers []float64, rawdata [][]float64) *Data {
	d := new(Data)
	//I prefer to copy the slice to avoid somebody changing it from outside
	d.dividers = make([]float64, len(dividers))
	for i, v := range dividers {
		d.dividers[i] = v
	}
	d.histo = make([]float64, len(dividers)-1)
	if rawdata != nil {
		d.total = len(rawdata)
		d.ReHisto(d.dividers, rawdata)
	}
	return d

}

//Adds the given data point(s) to the histogram
func (M *Data) AddData(point ...float64) {
	if D.normalized {
		D.UnNormalize()
	}
	for _, v := range point {
		for j, w := range M.dividers {
			if j == len(M.dividers)-1 {
				histo[j]++
				break
			}
			if w <= v && v < M.dividers[j+1] {
				histo[j]++
			}
		}
	}
	M.total += len(point)
	//if it was normalized, we should return it to that state
	if D.normalized {
		D.Normalize()
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
		return nil
	}
	n := D.total
	D.normalized = false
	if normalize {
		n = 1 / D.total
		D.normalized = true
	}

	D.histo = floats.Scale(n, D.histo)

}

//Copies the dividers of the histogram
func CopyDividers(dest ...[]float64) []float64 {
	d := getCopySlice(len(D.dividers), dest...)
	return floats.ScaleTo(d, 0, D.dividers)
}

func (D *Data) Copy(dest ...[]float64) []float64 {
	d := getCopySlice(len(D.histo), dest...)
	return floats.ScaleTo(d, 0, D.Histo)
}

func (D *Data) View() []float64 {
	return D.histo
}

func (D *Data) Add(a, b *Data) []float64 {
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
		D.histo[v] = a.histo[v] + b.histo[v]
	}
}

func (D *Data) Sub(a, b *Data, abs ...bool) []float64 {
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

		D.histo[v] = f(a.histo[v] - b.histo[v])
	}
}

func (D *Data) Sum() float64 {
	return floats.Sum(D.histo)
}

func (D *Data) ReHisto(dividers, rawdata []float64) {
	D.histo = stat.Histogram(dividers, rawdata)
}

func getCopySlice(N float64, dest ...[]float64) []float64 {
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
