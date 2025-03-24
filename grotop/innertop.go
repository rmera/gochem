package gro

import (
	"math"

	chem "github.com/rmera/gochem"
)

func sigmaepsilon2c6c2(sigma, epsilon float64) (c6 float64, c12 float64) {
	return 4 * sigma * math.Pow(epsilon, 6), sigma * 4 * (math.Pow(epsilon, 12))
}

type FF struct {
	Mol         *chem.Topology
	Bonds       []*Term
	Constraints []*Term
	Angles      []*Term
	Impropers   []*Term
	Dihedrals   []*Term
	Types       []*AtomType
	LJ          []*LJPair
}

type AtomType struct {
	Sigma   float64
	Epsilon float64
	C6      float64
	C12     float64
	AtNum   int
	Mass    float64
	Charge  float64
	Ptype   string
}

type LJPair struct {
	IDs     []int
	Names   []string
	Sigma   float64
	Epsilon float64
	C6      float64
	C12     float64
}

type Term struct {
	Functype   uint
	IDs        []int
	OneBased   int //0 if the indexes in Atom are 0-based, 1 if it's 1-based.
	K          float64
	Eq         float64
	Constraint bool
	Vsite      bool
	RB         []float64
}
