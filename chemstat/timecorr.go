package chemstat

import (
	"fmt"
	"math"
	"math/cmplx"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	"gonum.org/v1/gonum/dsp/fourier"
	"gonum.org/v1/gonum/stat"
)

func findClosestReal(list []complex128, index int) int {
	var min float64
	var ret int
	comp := real(list[index])
	min = math.Abs(comp - real(list[0]))
	for i, v := range list {
		if i == index {
			continue
		}
		if math.Abs(comp-real(v)) < min {
			min = math.Abs(comp - real(v))
			ret = i
		}
	}
	return ret

}

func cmplxMul(dst, b []complex128) {
	if len(dst) != len(b) {
		panic(fmt.Sprintf("complex multiplication of slices: Both slices should have the same len %d, %d", len(dst), len(b)))
	}
	for i, v := range b {
		dst[i] *= v
	}
}

func cmplxMulConj(dst, b []complex128) {
	if len(dst) != len(b) {
		panic(fmt.Sprintf("complex conjugate multiplication of slices: Both slices should have the same len %d, %d", len(dst), len(b)))
	}
	for i, v := range b {
		dst[i] *= cmplx.Conj(v)
	}
}

func cmplxRealScale(dst []complex128, sc float64) []complex128 {
	for i, v := range dst {
		//	fmt.Println(v, v*complex(sc, 1), sc) /////////////
		dst[i] = v * complex(sc, sc)
	}
	return dst
}

func reals(b []complex128, dst []float64) []float64 {
	for i, v := range b {
		dst[i] = real(v)
	}
	return dst
}

func conj(dst []complex128) []complex128 {
	for i, v := range dst {
		//fmt.Println("should be conjugates", dst[i], v, cmplx.Conj(v)) ///////////////////////////
		dst[i] = cmplx.Conj(v)
	}
	return dst
}

func mapfunc(t chem.Traj, coord *v3.Matrix, A chem.Atomer, c []float64, f func(c *v3.Matrix, A chem.Atomer) float64) []float64 {
	for err := t.Next(coord); ; err = t.Next(coord) {
		if err != nil {
			if _, ok := err.(chem.LastFrameError); ok {
				break
			}
			panic(err.Error())
		}
		c = append(c, f(coord, A))
	}
	return c
}

// Finds the cros-correlation function for the functions produced by f1 and f2 on trajectories t1 and t2
// which correspond to topologies A1 and A2, respectively. If all X1 and X2 (functions, topologies and trajs) are the same, you obtain the
// autocorrelation function. In that case, you can slightly lessen the computational cost by setting auto to true.
func MDCorrelation(t1, t2 chem.Traj, A1, A2 chem.Atomer, f1, f2 func(c *v3.Matrix, A chem.Atomer) float64, trjlen int, auto ...bool) []float64 {
	if trjlen == 0 {
		trjlen = 1e6
	}
	c1 := make([]float64, 0, trjlen)
	coord := v3.Zeros(t1.Len())
	c1 = mapfunc(t1, coord, A1, c1, f1)

	//fmt.Println("INITIAL DATA", len(c1), c1) /////////////
	c2 := make([]float64, len(c1))
	if len(auto) > 0 && auto[0] {
		copy(c2, c1)
	} else {
		coord = v3.Zeros(t1.Len())
		c2 = mapfunc(t2, coord, A2, c2, f2)
	}
	c1pad := make([]complex128, len(c1))
	c2pad := make([]complex128, len(c2))
	ret := make([]float64, 0, len(c1pad))
	ret = CrossCorrMem(c1, c2, c1pad, c2pad, ret) //, ret)
	//fmt.Println(len(c1), len(ret), maxlag)       //////////////////////
	return ret

}

// The difference between len(c1pad) and len(c1) (the pairs c1,c2conj and c1pad,c2conjpad must be of equal length)
// defines the lags to be used
func CrossCorrMem(c1, c2 []float64, c1pad, c2pad []complex128, dst ...[]float64) []float64 {

	//	lags := len(c1pad) - len(c1)
	var ret []float64
	if len(dst) == 0 || len(dst[0]) > 0 { //if you give a slice, you can the cap, but len must be 0
		ret = make([]float64, 0, len(c1pad))
	} else {
		ret = dst[0]
	}

	//we start preparing the data
	c1mean := stat.Mean(c1, nil)
	c2mean := stat.Mean(c2, nil)
	c1std := stat.StdDev(c1, nil)
	c2std := stat.StdDev(c2, nil)
	if len(c1pad) != 2*len(c1) {
		c1pad = make([]complex128, 2*len(c1))
	}
	if len(c2pad) != 2*len(c2) {
		c2pad = make([]complex128, 2*len(c2))
	}
	for i, v := range c1 {
		c1pad[i] = complex(v-c1mean, 0)
		c2pad[i] = complex(c2[i]-c2mean, 0)
	}
	//	println("partimos")
	f := fourier.NewCmplxFFT(len(c1pad)) //(len(c1pad))
	f.Coefficients(c1pad, c1pad)         //(c1pad, c1pad)
	f.Coefficients(c2pad, c2pad)         //(c2pad, c2pad)
	cmplxMulConj(c1pad, c2pad)
	f.Sequence(c1pad, c1pad)

	cmplxRealScale(c1pad, (1.0 / float64(len(c1pad)))) //normalization of the FFT

	//c1pad = c1pad[0 : len(c1pad)-1]
	center := len(c1pad) / 2

	for _, v := range c1pad[:center] {
		ret = append(ret, real(v))
	}
	for _, v := range c1pad[center:] {
		ret = append(ret, (real(v))) // / float64(i)))
	}
	//	log.Println("LEN RET", len(ret)) ////////////////////
	for i, v := range ret {
		ret[i] = v / (c1std * c2std) / float64(len(c1))
	}
	//	fmt.Println(ret[:100])         ////////////
	//	fmt.Println(ret[len(ret)-38:]) //////////
	return ret
}

func RMSDCorrFunc(refcoord *v3.Matrix, indexes ...[]int) func(c *v3.Matrix, A chem.Atomer) float64 {
	return func(c *v3.Matrix, A chem.Atomer) float64 {
		var f float64
		var err error
		if len(indexes) > 0 {
			f, err = chem.RMSD(c, refcoord, indexes[0], indexes[0])
		} else {
			f, err = chem.RMSD(c, refcoord)
		}
		if err != nil {
			panic(err.Error())
		}
		return f
	}
}
