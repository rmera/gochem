package protcap

import (
	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

type ComplexBBCap struct {
	cp *capper
}

func NewComplexBBCap() *ComplexBBCap {
	r := new(ComplexBBCap)
	r.cp = newCapper()
	return r
}

func (CC *ComplexBBCap) Cap(c *v3.Matrix, mol *chem.Topology, tocapN, tocapC int, putLast bool) (*v3.Matrix, *chem.Topology) {
	return BBHCap(c, mol, tocapN, tocapC, putLast, CC.cp)

}

// CapBackBoneNoRef takes a sub-peptide, represented by c and mol, of a larger peptide and caps the N- side of the first residue
// with MolID tocapN and the C- side of the residue with MolID tocapC (they can both be the same residue). If you only want to cap
// one of both sides, you can give -1 for tocapN or tocapC to avoid adding an N-cap or C-cap, respectively.
// By default, the capping atom will be added to the end of the corresponding residue, but you can put it at the end of the
// c and mol by giving a true in putLast.
func BBHCap(c *v3.Matrix, mol *chem.Topology, tocapN, tocapC int, putLast bool, cp ...*capper) (*v3.Matrix, *chem.Topology) {
	// new method to cap C=O & C-N. For C=O, take the vectors C-O & C-CA, make a copy of the C-CA one. Take the cross product between C-O and C-CA, and rotate the copy around the cross product until it is in the right position. (so ~60 degrees, I'd think).
	var iN, iC, iCAN, iCAO, iO, iH int
	Zero := v3.Zeros(1)
	V1 := v3.Zeros(1)
	V2 := v3.Zeros(1)
	Cross := v3.Zeros(1)
	Ncap := v3.Zeros(1)
	Ccap := v3.Zeros(1)

	for i := 0; i < mol.Len(); i++ {
		a := mol.Atom(i)
		if tocapN == a.MolID { //slices.Contains(tocapN, a.MolID) {
			switch a.Name {
			case "N":
				iN = i
			case "H":
				iH = i
			case "CA":
				iCAN = i
			case "CD":
				if a.MolName == "PRO" {
					iH = i //proline doesn't have NH since the CD carbon form the lateral chain is
					//bonded to the N.
				}
			}
		}
		if tocapC == a.MolID { //slices.Contains(tocapC, a.MolID) {
			switch a.Name {
			case "C":
				//	println("C!", i)
				iC = i
			case "O":
				iO = i
			case "CA":
				iCAO = i

			}
		}
	}
	rotangle := 120 * chem.Deg2Rad
	//The N-Capping
	placekeep := -1
	if tocapN > 0 {
		nref := mol.Atom(iCAN)
		place := nref.MolID
		added := 0
		if putLast {
			place = -1 * place
		} else {
			placekeep = place
		}

		if len(cp) == 0 {

			Ncap = capGeo(c, Ncap, Zero, V1, V2, Cross, iN, iH, iCAN, -1*rotangle)
			nc := new(chem.Atom)
			nc.Name = "HNA"
			nc.Symbol = "H"
			nc.ID = mol.Len()
			nc.SetIndex(nc.ID - 1)
			nc.MolName = nref.MolName
			nc.MolID = nref.MolID
			nc.Chain = nref.Chain
			c, mol = insertInMol(Ncap, c, []*chem.Atom{nc}, mol, place)
			added = 1
		} else {
			ats, tc := cp[0].NCap(c, iN, iH, iCAN, nref)

			c, mol = insertInMol(tc, c, ats.Atoms, mol, place)
			added = 3
		}
		//must update the C-cap indexes to consider the atoms we added in the middle
		//as N-cap. placekeep will only be >0 if we add N-cap and we add it to the middle.
		if placekeep > 0 {
			if mol.Atom(iC).MolID > placekeep { //this check shouldn't be needed, ncap should always
				//come before cap so the condition should always be true.
				iC += added
				iCAO += added
				iO += added
			}
		}

	}

	//The C-Capping
	if tocapC > 0 {
		cref := mol.Atom(iCAO)
		place := cref.MolID
		if putLast {
			place = -1 * place
		}

		if len(cp) == 0 {
			Ccap = capGeo(c, Ccap, Zero, V1, V2, Cross, iC, iO, iCAO, -1*rotangle)
			//we now create the capping atoms to add to the topology
			cc := new(chem.Atom)
			cc.Name = "HCP"
			cc.Symbol = "H"
			cc.ID = mol.Len() + 1
			cc.SetIndex(cc.ID - 1)
			cc.MolID = cref.MolID
			cc.MolName = cref.MolName
			cc.Chain = cref.Chain
			c, mol = insertInMol(Ccap, c, []*chem.Atom{cc}, mol, place)
		} else {

			ats, tc := cp[0].CCap(c, iC, iO, iCAO, cref)
			c, mol = insertInMol(tc, c, ats.Atoms, mol, place)
		}

	}

	return c, mol
}

// inserts in in c and ina in mol as the last element of the
// residue with ID molID. If molID is not given or absent from mol, it will insert the atom at the end
func insertInMol(in, c *v3.Matrix, inat []*chem.Atom, mol *chem.Topology, molID int) (*v3.Matrix, *chem.Topology) {
	posbefore := mol.Len() - 1
	//	first := mol.Atom(0).MolID //this is to ensure that if the given molID is the first residue
	//we don't pick the first atom in that residue, but the last
	for i := 0; i < mol.Len(); i++ {
		a := mol.Atom(i)
		if a.MolID == molID+1 { // && a.MolID != first {
			posbefore = i - 1 //here we are already in the residue we want.
			//posbefore is the position before that, so we subtract one
			break
		}
	}

	pre := c.View(0, 0, posbefore+1, 3) //posbefore is the index (0-based) of the last element, so posbefore+1 is the
	//number of elements to take.
	ret := v3.Zeros(c.Len() + in.Len())
	ret.Stack(pre, in)
	start := posbefore + in.Len() + 1
	if posbefore < c.Len()-1 {
		post := c.View(posbefore+1, 0, c.Len()-posbefore-1, 3)
		for i := start; i < ret.Len(); i++ {
			x := post.At(i-start, 0)
			y := post.At(i-start, 1)
			z := post.At(i-start, 2)
			ret.Set(i, 0, x)
			ret.Set(i, 1, y)
			ret.Set(i, 2, z)
		}
	}
	//now the atoms
	ats := make([]*chem.Atom, 0, mol.Len()+len(inat))
	ats = append(ats, mol.Atoms[0:posbefore+1]...)
	//	inat2 := make([]*chem.Atom, 0, len(inat))
	for i, v := range inat {
		w := new(chem.Atom)
		w.Copy(v)
		w.ID = posbefore + 1 + i
		w.SetIndex(posbefore + i)
		w.MolID = molID
		if molID < 0 {
			w.MolID = -1 * molID
		}
		ats = append(ats, w)
	}
	//	ats = append(ats, inat2...)
	if posbefore < c.Len()-1 {
		ats = append(ats, mol.Atoms[posbefore+1:]...)
	}
	tret := chem.NewTopology(0, 1, ats)
	return ret, tret
}

func capGeo(c, Cap, Zero, V1, V2, Cross *v3.Matrix, izero, iv1, iv2 int, angle float64) *v3.Matrix {
	Zero.Copy(c.VecView(izero))
	c.SubVec(c, Zero)
	Cap.Copy(c.VecView(iv1))
	V1.Copy(Cap)
	V2.Copy(c.VecView(iv2))
	Cross.Cross(V1, V2)
	Cap = chem.Rotate(V1, Cap, Cross, angle)
	Cap.AddVec(Cap, Zero)
	c.AddVec(c, Zero)
	return Cap
}
