package clash

import (
	"fmt"
	"log"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/top"
	v3 "github.com/rmera/gochem/v3"

	"gonum.org/v1/gonum/floats"
)

func bondRotate(coord *v3.Matrix, at1, at2 int, angle float64, torotate []int, temp ...*v3.Matrix) (*v3.Matrix, error) {
	var tmp *v3.Matrix
	if len(temp) > 0 && temp[0].Len() == len(torotate) {
		tmp = temp[0]
	} else {
		tmp = v3.Zeros(len(torotate))
	}
	a1 := coord.VecView(at1)
	a2 := coord.VecView(at2)
	tmp.SomeVecs(coord, torotate)
	nc, err := chem.RotateAbout(tmp, a1, a2, angle)
	if err != nil {
		return nil, err
	}
	coord.SetVecs(nc, torotate)
	return coord, nil
}

func GeometryOnlyHighestOverlap(testtop, clashtop chem.Atomer, FF top.FF) func(*v3.Matrix, *v3.Matrix) (float64, bool, [2]int) {
	return func(t *v3.Matrix, c *v3.Matrix) (float64, bool, [2]int) {
		return HighestOverlap(t, c, testtop, clashtop, FF)
	}

}

func pairradiioverlap(distance float64, t, c int, ttop, ctop chem.Atomer, FF top.FF) (float64, error) {
	tname := ttop.Atom(t).Name
	cname := ctop.Atom(c).Name
	radsum := FF.VdWForPair(tname, cname)
	if radsum <= 0 {
		return 1 / distance, fmt.Errorf("Couldn't find radii for atoms %s-%s", tname, cname)
	}
	return radsum - distance, nil
}

func HighestOverlap(test, clash *v3.Matrix, testtop, clashtop chem.Atomer, FF top.FF) (over float64, in bool, indexes [2]int) {
	overIn := 0.0
	dvec := v3.Zeros(1)
	dt := 0.0
	var a1, a2 *v3.Matrix
	var iindexes [2]int
	for i := 0; i < test.Len(); i++ {
		for j := 0; j < i; j++ {
			a1 = test.VecView(i)
			a2 = test.VecView(j)
			dvec.SubVec(a1, a2)
			dt = dvec.Norm(2)
			ov, err := pairradiioverlap(dt, i, j, testtop, testtop, FF)
			if err != nil {
				log.Println(err)
				continue
			}
			if ov > overIn {
				overIn = dt
				iindexes[0] = i
				iindexes[1] = j
			}
		}
	}
	for i := 0; i < test.Len(); i++ {
		for j := 0; j < clash.Len(); j++ {
			a1 = test.VecView(i)
			a2 = clash.VecView(j)
			dvec.SubVec(a1, a2)
			dt = dvec.Norm(2)
			ov, err := pairradiioverlap(dt, i, j, testtop, clashtop, FF)
			if err != nil {
				log.Println(err)
				continue
			}

			if ov < over {
				over = ov
				indexes[0] = i
				indexes[1] = j
			}
		}
	}
	if overIn < over {
		in = true
		indexes = iindexes
		over = overIn
	}
	return
}

func LowestDist(test, clash *v3.Matrix) (dist float64, in bool, indexes [2]int) {
	dist = 999999999999999
	distIn := 999999999999.0
	dvec := v3.Zeros(1)
	dt := 0.0
	var a1, a2 *v3.Matrix
	var iindexes [2]int
	for i := 0; i < test.Len(); i++ {
		for j := 0; j < i; j++ {
			a1 = test.VecView(i)
			a2 = test.VecView(j)
			dvec.SubVec(a1, a2)
			dt = dvec.Norm(2)
			if dt < distIn {
				distIn = dt
				iindexes[0] = i
				iindexes[1] = j
			}

		}

	}
	for i := 0; i < test.Len(); i++ {
		for j := 0; j < clash.Len(); j++ {
			a1 = test.VecView(i)
			a2 = clash.VecView(j)
			dvec.SubVec(a1, a2)
			dt = dvec.Norm(2)
			if dt < dist {
				dist = dt
				indexes[0] = i
				indexes[1] = j
			}

		}
	}
	if distIn < dist {
		in = true
		indexes = iindexes
		dist = distIn
	}
	return
}

// a q&d brute-force central-difference
func AngleFuncGrad(test, clash *v3.Matrix, axes, rotated [][]int, f func(*v3.Matrix, *v3.Matrix) (float64, bool, [2]int), epsilon ...float64) ([]float64, error) {
	var e float64 = 2 * chem.Deg2Rad
	if len(epsilon) > 0 {
		e = epsilon[0]
	}
	clone := v3.Zeros(test.Len())
	grad := make([]float64, 0, len(axes))
	//	dist,_,_:=LowestDist(test,clash)
	for i, v := range axes {
		clone.Copy(test)
		rot, err := bondRotate(clone, v[0], v[1], e, rotated[i])
		if err != nil {
			return nil, err
		}
		dpos, _, _ := f(rot, clash)
		clone.Copy(test)
		rot, err = bondRotate(clone, v[0], v[1], -1*e, rotated[i])
		if err != nil {
			return nil, err
		}
		dneg, _, _ := f(rot, clash)
		g := (dpos - dneg) / 2 * e //central difference
		grad = append(grad, g)
	}
	floats.Scale(floats.Sum(grad), grad) //normalization
	return grad, nil
}

// A q&d brute-force step in the angles space (the angles defined by axes and rotated) with step step[0], or
// 10 degrees if not given, in the direction of the gradient grad, which is assumed to be normalized.
// If you want to go in the opposite direction, simply give a negative step.
func AngleStep(test *v3.Matrix, axes, rotated [][]int, grad []float64, step ...float64) (*v3.Matrix, error) {
	var s float64 = 10 * chem.Deg2Rad
	var err error
	if len(step) > 0 {
		s = step[0]
	}
	clone := v3.Zeros(test.Len())
	clone.Copy(test)
	//	dist,_,_:=LowestDist(test,clash)
	for i, v := range axes {
		clone, err = bondRotate(clone, v[0], v[1], s*grad[0], rotated[i])
		if err != nil {
			return nil, err
		}
	}
	return clone, nil
}
