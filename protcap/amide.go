package protcap

import (
	"strings"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

// I think this is small enough and general enough to be worth hardcoding. It's a formamide molecule, optimized at the
// GFN2 level.
const capstr string = "6\n\nC           -4.92794702449458        1.89045342665579       -0.07921100123582\nO           -4.70801763850576        3.04272114512478       -0.33847635299908\nN           -4.42554476863113        1.23149344587307        0.98228022653061\nH           -4.63126053618251        0.26640246053375        1.16518717959352\nH           -3.81743860256653        1.72542300625019        1.61607772779290\nH           -5.57786142961950        1.24519651556242       -0.69948777968212\n"

const CN float64 = 1.35 //C-N distance in amide (aprox) in A

type capper struct {
	mol       *chem.Molecule
	coord     *v3.Matrix
	coord2    *v3.Matrix
	tmp       *v3.Matrix
	tmp2      *v3.Matrix
	tmpNCap   *v3.Matrix
	tmpCCap   *v3.Matrix
	HNReplace int
	HCReplace int
	N         int
	C         int
	O         int
	HN        int //the no-replaced one.
}

func newCapper() *capper {
	r := new(capper)
	mol, err := chem.XYZRead(strings.NewReader(capstr)) //as this is hardcoded, thre really shouldn't be errors.
	//I can just panic if that happens as it is surely a bug.
	if err != nil {
		panic(err)
	}
	r.tmp = v3.Zeros(1)
	r.tmp2 = v3.Zeros(1)
	r.tmpCCap = v3.Zeros(3)
	r.tmpNCap = v3.Zeros(3)
	r.mol = mol
	r.coord = mol.Coords[0]
	r.coord2 = v3.Zeros(r.coord.Len())
	r.HNReplace = 4
	r.HCReplace = 5
	r.C = 0
	r.N = 2
	r.O = 1
	r.HN = 3
	return r
}

// Returns a matrix and topology containing the atoms and proper positions to C-cap the residue with coordinates
// in co, with indexes for the atoms C,O and CA given. firstID is the ID wanted for the first
// atom in the capping. sample is an atom from the residue to be capped, from which to get things like
// MolID, MolName and Chain.
func (c *capper) CCap(co *v3.Matrix, C, O, CA int, sample *chem.Atom) (*chem.Topology, *v3.Matrix) {
	c.coord2.Copy(c.coord)

	//c.cap will superimpose atoms in such a way that we need to scale our C-H  bond to be as long as an C-N bond
	//so the superposition worke better. The displaced H atom will not be part of the returned molecule anyway.

	c.StretchToDist(c.HCReplace, c.C, CN, []int{c.HCReplace}, c.tmp2)
	return c.cap(co, c.tmpCCap, []int{c.N, c.HN, c.HNReplace}, []int{c.HCReplace, c.C, c.O}, []int{CA, C, O}, sample)
}

// Returns a matrix and topology containing the atoms and proper positions to N-cap the residue with coordinates
// in co, and with indexes for the atoms N,HN and CA given. firstID is the ID wanted for the first
// atom in the capping. sample is an atom from the residue to be capped, from which to get things like
// MolID, MolName and Chain.

func (c *capper) NCap(co *v3.Matrix, N, HN, CA int, sample *chem.Atom) (*chem.Topology, *v3.Matrix) {
	c.coord2.Copy(c.coord)
	//c.cap will superimpose atoms in such a way that we need to scale our N-H bond to be as long as an C-N bond
	//so the superposition worke better. The displaced H atom will not be part of the returned molecule anyway.
	c.StretchToDist(c.HNReplace, c.N, CN, []int{c.HNReplace}, c.tmp2)
	return c.cap(co, c.tmpNCap, []int{c.C, c.O, c.HCReplace}, []int{c.N, c.HN, c.HNReplace}, []int{N, HN, CA}, sample)

}

// This is the more general function used by CCap and NCap.
func (c *capper) cap(co *v3.Matrix, tmp *v3.Matrix, capindexes, supindexes, resindexes []int, sample *chem.Atom) (*chem.Topology, *v3.Matrix) {
	c.coord2.Copy(c.coord)
	//	c.StretchToDist(c.HNReplace, c.N, CN, []int{c.HNReplace}, c.tmp2)
	var err error
	c.coord2, err = chem.Super(c.coord2, co, supindexes, resindexes)
	if err != nil {
		panic("Superposition failed for capping:" + err.Error()) //NOTE: Consider changing to return an error
	}
	ats := []*chem.Atom{c.mol.Atom(capindexes[0]), c.mol.Atom(capindexes[1]), c.mol.Atom(capindexes[2])}
	//	t := v3.Zeros(len(capindexes))
	tmp.SomeVecs(c.coord2, capindexes)
	for _, v := range ats {
		//	v.ID = firstID + i
		v.SetIndex(v.ID - 1)
		v.MolID = sample.MolID
		v.MolName = sample.MolName
		v.Name = v.Name + "C"
		v.Chain = sample.Chain
		if len(v.Name) > 3 {
			v.Name = v.Name[:3] //not to keep adding 'C's very time this is called.
		}

	}
	return chem.NewTopology(0, 1, ats), tmp
}

// Moves all atoms in coord2 with indexes present in topull along the axis between the atoms
// i and j in coord2, until the distance between i and j (if one of them is include in topull) is
// dist. tmp is given for temp storage.
func (c *capper) StretchToDist(i, j int, dist float64, topull []int, tmp *v3.Matrix) float64 {
	c.tmp.Sub(c.coord2.VecView(j), c.coord2.VecView(i))
	tdist := c.tmp.Norm(2)
	c.tmp.Unit(c.tmp)
	newmag := dist - tdist
	c.tmp.Scale(newmag, c.tmp)
	tmp.SomeVecs(c.coord2, topull)
	tmp.AddVec(tmp, c.tmp)
	c.coord2.SetVecs(tmp, topull)
	return tdist
}

// Returns a copy of c.
func (c *capper) Copy() *capper {
	r := newCapper()
	r.coord.Copy(r.coord)
	r.coord2.Copy(r.coord2)
	return r

}
