package histo

import (
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/stat"
)

type Data struct {
	doHisto   bool
	rawdata   []float64
	dividers  []float64
	histogram []float64
}

func (D *Data) AddData(point ...float64) {
	D.rawdata = append(D.rawdata, point...)
	if D.doHisto {
		D.histogram = stat.Histogram(D.rawdata, D.dividers)
	}
}

type MData struct {
	normalized bool
	total      int
	dividers   []float64
	histo      []float64
}

func (M *MData) AddData(point ...float64) {
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

func (D *MData) Normalized() bool {
	return D.normalized
}

func (D *MData) Normalize() {
	if D.total == 0 {
		return nil
	}
	//we don't waste time if the histogram is already normalized.
	if !D.normalized {
		D.histo = floats.Scale(1/D.total, D.histo)
		D.normalized = true
	}
}

func (D *MData) UnNormalize() {
	if D.total == 0 {
		return nil
	}
	if D.normalized {
		D.histo = floats.Scale(D.total, D.histo)
		D.normalized = false
	}
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

func CopyDividers() []float64 {
	d := getCopySlice(len(D.dividers), dest...)
	return floats.ScaleTo(d, 0, D.dividers)

}

func (D *MData) Copy(dest ...[]float64) []float64 {
	d := getCopySlice(len(D.histo), dest...)
	return floats.ScaleTo(d, 0, D.Histo)
}

func (D *MData) View() []float64 {
	return D.histo
}

func (D *MData) ReHisto(dividers, rawdata []float64) {
	D.histo = stat.Histogram(dividers, rawdata)
}
