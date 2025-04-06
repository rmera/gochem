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

// returns a new and empty (but with some values set to defaults)
// FF object. SigmaEpsilon is true if LJ terms are expressed as sigma/epsilon,
// false for C6/C12
func NewFF(mol *chem.Topology, SigmaEpsilon ...bool) *FF {
	se := true //sigma-epsilon true is the form used by Martini3.
	if len(SigmaEpsilon) > 0 {
		se = SigmaEpsilon[0]
	}
	ret := new(FF)
	ret.Mol = mol
	ret.SigmaEpsilon = se
	ret.currentHeader = "NOHEADER"
	return ret
}

// Returns the nummber of atoms in the topology.
func (F *FF) Len() int {
	return F.Mol.Len()
}

// Returns a copy of the receiver
func (F *FF) Copy() *FF {
	F2 := new(FF)
	mol := chem.NewTopology(0, 1)
	for i := 0; i < F.Mol.Len(); i++ {
		a := new(chem.Atom)
		a.Copy(F.Mol.Atom(i))
		mol.AppendAtom(a)
	}
	F2.Mol = mol
	F2.SigmaEpsilon = F.SigmaEpsilon
	F2.currentHeader = F.currentHeader
	F2.Bonds = copyTerms(F.Bonds)
	F2.Constraints = copyTerms(F.Constraints)
	F2.Angles = copyTerms(F.Angles)
	F2.Impropers = copyTerms(F.Impropers)
	F2.Dihedrals = copyTerms(F.Dihedrals)
	F2.ATypes = make([]*AtomType, 0, len(F.ATypes))
	for _, v := range F.ATypes {
		a := new(AtomType)
		a.Copy(v)
		F2.ATypes = append(F2.ATypes, a)
	}
	F2.LJ = make([]*LJPair, 0, len(F.LJ))
	for _, v := range F.LJ {
		a := new(LJPair)
		a.Copy(v)
		F2.LJ = append(F2.LJ, a)
	}
	F2.Exclusions = make([][]int, 0, len(F.Exclusions))
	for _, v := range F.Exclusions {
		a := make([]int, len(v))
		copy(a, v)
		F2.Exclusions = append(F2.Exclusions, a)
	}
	F2.VSites = make([]*VSite, 0, len(F.VSites))
	for _, v := range F.VSites {
		a := new(VSite)
		a.Copy(v)
		F2.VSites = append(F2.VSites, a)
	}
	return F2
}

func copyTerms(T []*Term) []*Term {
	ret := make([]*Term, 0, len(T))
	for _, v := range T {
		t2 := new(Term)
		t2.Copy(v)
		ret = append(ret, t2)
	}
	return ret
}

type AtomType struct {
	SigmaEpsilon bool
	Name         string
	C6           float64
	C12          float64
	AtNum        int
	Mass         float64
	Charge       float64
	Ptype        string
}

// Copy puts a copy of B in the receiver
func (A *AtomType) Copy(B *AtomType) {
	A.SigmaEpsilon = B.SigmaEpsilon
	A.Name = B.Name
	A.C6 = B.C6
	A.C12 = B.C12
	A.AtNum = B.AtNum
	A.Mass = B.Mass
	A.Charge = B.Charge
	A.Ptype = B.Ptype
}

func (A *AtomType) equal(B any) bool {
	b := B.(*AtomType)
	return A.Name == b.Name
}

type LJPair struct {
	//	IDs   []int
	SigmaEpsilon bool
	Names        []string
	FuncType     int
	C6           float64
	C12          float64
}

// Copies B into the receiver
func (A *LJPair) Copy(B *LJPair) {
	A.SigmaEpsilon = B.SigmaEpsilon
	if len(A.Names) != len(B.Names) {
		A.Names = make([]string, len(B.Names))
	}
	copy(A.Names, B.Names)
	A.C6 = B.C6
	A.C12 = B.C12
	A.FuncType = B.FuncType
}

func (A *LJPair) equal(B any) bool {
	b := B.(*LJPair)
	sym := A.Names[0] == b.Names[0] && A.Names[1] == b.Names[1]
	asym := A.Names[1] == b.Names[0] && A.Names[0] == b.Names[1]
	return sym || asym
}

type VSite struct {
	ID       int
	N        int //0 for virtual_sistesn
	FuncType int
	Atoms    []int
	Factors  []float64
}

// Copies B into the receiver
func (A *VSite) Copy(B *VSite) {
	A.ID = B.ID
	A.N = B.N
	A.FuncType = B.FuncType
	if len(A.Atoms) != len(B.Atoms) {
		A.Atoms = make([]int, len(B.Atoms))
	}
	copy(A.Atoms, B.Atoms)
	if len(A.Factors) != len(B.Factors) {
		A.Factors = make([]float64, len(B.Factors))
	}
	copy(A.Factors, B.Factors)

}

type Term struct {
	FuncType   uint
	IDs        []int
	OneBased   int //1 if the indexes in Atom are 0-based, 0 if it's 1-based.
	K          float64
	Eq         float64
	Constraint bool
	Vsite      bool
	RB         []float64
}

func (A *Term) Copy(B *Term) {
	if len(A.IDs) != len(B.IDs) {
		A.IDs = make([]int, len(B.IDs))
	}
	copy(A.IDs, B.IDs)
	A.OneBased = B.OneBased
	A.FuncType = B.FuncType
	A.K = B.K
	A.Eq = B.Eq
	A.Constraint = B.Constraint
	A.Vsite = B.Vsite
	if len(A.RB) != len(B.RB) {
		A.RB = make([]float64, len(B.RB))
	}
	copy(A.RB, B.RB)
}
