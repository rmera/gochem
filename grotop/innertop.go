package gro

import (
	"math"

	chem "github.com/rmera/gochem"
)

func sigmaepsilonToc6c2(sigma, e float64) (c6 float64, c12 float64) {
	return 4 * e * math.Pow(sigma, 6), e * 4 * (math.Pow(sigma, 12))
}

func c6c12ToSigmaepsilon(c6, c12 float64) (sigma float64, epsilon float64) {
	return math.Pow(c6/c12, (1 / 6)), math.Pow(c6, 2) / (4 * c12)
}

type FF struct {
	SigmaEpsilon  bool //are LJ terms using sigma/epsilon, or C6/C12?
	currentHeader string
	Mol           *chem.Topology
	Bonds         []*Term
	Constraints   []*Term
	Angles        []*Term
	Impropers     []*Term
	Dihedrals     []*Term
	VSites        []*VSite
	ATypes        []*AtomType
	LJ            []*LJPair
	Exclusions    [][]int
}

func NewFF(mol *chem.Topology, SigmaEpsilon ...bool) *FF {
	se := true //sigma-epsilon true is the form used by Martini3.
	if len(SigmaEpsilon) > 0 {
		se = SigmaEpsilon[0]
	}
	ret := new(FF)
	ret.Mol = mol
	ret.SigmaEpsilon = se
	return ret
}

type AtomType struct {
	Name   string
	C6     float64
	C12    float64
	AtNum  int
	Mass   float64
	Charge float64
	Ptype  string
}

type LJPair struct {
	//	IDs   []int
	Names    []string
	FuncType int
	C6       float64
	C12      float64
}

type VSite struct {
	ID       int
	N        int //0 for virtual_sistesn
	FuncType int
	Atoms    []int
	Factors  []float64
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
