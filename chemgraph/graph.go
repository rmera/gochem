package chemgraph

import (
	"fmt"
	"math"

	chem "github.com/rmera/gochem"
	"gonum.org/v1/gonum/graph"
)

type Atom struct {
	*chem.Atom
	Bonds  []*Bond
	IDFunc func(*Atom) int64
}

func (A *Atom) ID() int64 {
	if A.IDFunc == nil {
		return int64(A.Index())
	}
	return A.IDFunc(A)
}

func (A *Atom) AtID() int {
	return A.Atom.ID
}

type Bond struct {
	*chem.Bond
	At1, At2   *Atom
	Weightfunc func(*Bond) float64
}

func (B *Bond) Weight() float64 {
	if B.Weightfunc == nil {
		return B.Energy
	}
	return B.Weightfunc(B)
}

func (B *Bond) From() graph.Node {
	return B.At1
}

func (B *Bond) To() graph.Node {
	return B.At2
}

// I'll try to just switch them in place, as bonds are not directional
// look here if you have issues
func (B *Bond) ReversedEdge() graph.Edge {
	B.At1, B.At2 = B.At2, B.At1
	return B
}

type Bonds []*Bond

func (B Bonds) Len() int {
	return len(B)
}
func (B Bonds) Contains(index int) bool {
	for _, b := range B {
		if b.Index == index {
			return true
		}
	}
	return false
}

// Implements gonum.graph.Nodes
type Atoms struct {
	Atoms []*Atom
	curr  int
}

func (A *Atoms) Len() int {
	return len(A.Atoms)
}
func (A *Atoms) Reset() {
	A.curr = 0
}
func (A *Atoms) Next() bool {
	if A.curr >= len(A.Atoms) {
		return false
	}
	A.curr++
	return true
}
func (A *Atoms) Node() graph.Node {
	return A.Atoms[A.curr]
}

// implements Gonum graph.Graph and graph.Weighted interfaces
type Topology struct {
	*chem.Topology
	Bonds
	*Atoms
}

func (T *Topology) Nodes() *Atoms {
	if T.Atoms.Len() == 0 {
		panic("Topology has no atoms")
	}
	return T.Atoms
}

func (T *Topology) From(id int64) graph.Nodes {
	ret := make([]*Atom, 0)
	for _, b := range T.Bonds {
		///undirected graph
		if b.From().ID() == id || b.To().ID() == id {
			ret = append(ret, b.To().(*Atom))
		}
	}
	return &Atoms{curr: 0, Atoms: ret}
}

func (T *Topology) HasEdgeBetween(id1, id2 int64) bool {
	for _, b := range T.Bonds {
		if (b.From().ID() == id1 && b.To().ID() == id2) || (b.From().ID() == id2 && b.To().ID() == id1) {
			return true
		}
	}
	return false
}

func (T *Topology) WeightedEdgeBetween(id1, id2 int64) *Bond {
	for _, b := range T.Bonds {
		if (b.From().ID() == id1 && b.To().ID() == id2) || (b.From().ID() == id2 && b.To().ID() == id1) {
			return b
		}
	}
	return nil
}

func (T *Topology) Edge(id1, id2 int64) graph.Edge {
	for _, b := range T.Bonds {
		//I'm making the graph always undirected
		if (b.From().ID() == id1 && b.To().ID() == id2) || (b.From().ID() == id2 && b.To().ID() == id1) {
			return b
		}
	}
	return nil
}

func (T *Topology) WeightedEdge(id1, id2 int64) graph.Edge {
	return T.Edge(id1, id2)
}

func (T *Topology) Weight(id1, id2 int64) (w float64, ok bool) {
	if id1 == id2 {
		return 0.0, true
	}
	b := T.Edge(id1, id2)
	if b == nil {
		return -1, false
	}
	return b.(*Bond).Weight(), true
}

func atomID(Ats []*Atom, id int64) *Atom {
	for _, v := range Ats {
		if v.ID() == id {
			return v
		}
	}
	return nil
}

func bondContains(B Bonds, index int) bool {
	for _, b := range B {
		if b.Index == index {
			return true
		}
	}
	return false
}

// I might add an interafece to avoid using chem.Molecule
func TopologyFromChem(mol *chem.Molecule, IDFunc func(*Atom) int64, weightfunc func(*Bond) float64) *Topology {
	mol.FillIndexes()
	b := make([]*Bond, 0)
	a := make([]*Atom, 0)
	if IDFunc == nil {
		IDFunc = func(A *Atom) int64 { return int64(A.Index()) }
	}
	if weightfunc == nil {
		weightfunc = func(B *Bond) float64 { return 1 / math.Abs(B.Energy) }
	}
	for i := 0; i < mol.Len(); i++ {
		at := mol.Atom(i)
		a = append(a, &Atom{Atom: at, IDFunc: IDFunc})
		for _, v := range at.Bonds {
			if !bondContains(b, v.Index) {
				nb := &Bond{Bond: v, At1: &Atom{Atom: v.At1, IDFunc: IDFunc}, At2: &Atom{Atom: v.At2, IDFunc: funcID}, Weightfunc: weightfunc}
				b = append(b, nb)
			}
		}
	}

	for i, v := range b {
		At1 := atomID(a, v.At1.ID())
		At2 := atomID(a, v.At2.ID())
		if At1 == nil || At2 == nil {
			panic(fmt.Sprintf("TopologyFromChem: Bond %d has at least one non-existent atom", i))
		}
		At1.Bonds = append(At1.Bonds, v)
		At2.Bonds = append(At2.Bonds, v)
	}
	return &Topology{Topology: mol.Topology, Bonds: Bonds(b), Atoms: &Atoms{Atoms: a, curr: 0}}
}
