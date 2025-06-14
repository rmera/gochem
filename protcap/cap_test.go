package protcap

import (
	"fmt"
	"testing"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
)

const dir string = "../test"

func errql(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
	}
}

func TestCComplexCap(Te *testing.T) {
	mol, err := chem.PDBFileRead(dir + "/2c9vIOH.pdb")
	errql(Te, err)
	sq := [][]int{{13, 14, 15}, {36, 37, 38}}
	ch := []string{"A"}
	cp := NewComplexBBCap()
	for i, v := range sq {
		pep := chem.Molecules2Atoms(mol, v, ch)
		m1 := chem.NewTopology(0, 1)
		m1.SomeAtoms(mol, pep)
		c := v3.Zeros(len(pep))
		c.SomeVecs(mol.Coords[0], pep)
		ccoord, cmol := cp.Cap(c, m1, v[0], v[len(v)-1], true)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dPutLast.pdb", dir, i), ccoord, cmol, nil)
		cccoord, ccmol := cp.Cap(c, m1, v[0], v[len(v)-1], false)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dPutMed.pdb", dir, i), cccoord, ccmol, nil)
		ccccoord, cccmol := BBHCap(c, m1, v[0], v[len(v)-1], true)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dHPutLast.pdb", dir, i), ccccoord, cccmol, nil)
		cccccoord, ccccmol := BBHCap(c, m1, v[0], v[len(v)-1], false)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dHPutMed.pdb", dir, i), cccccoord, ccccmol, nil)

		wcord, wmol := cp.Cap(c, m1, -1, v[len(v)-1], true)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dOnlyC.pdb", dir, i), wcord, wmol, nil)
		wwcord, wwmol := cp.Cap(c, m1, v[0], -1, false)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dOnlyNPutMed.pdb", dir, i), wwcord, wwmol, nil)
		wwwcord, wwwmol := cp.Cap(c, m1, v[0], -1, true)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dOnlyNPutLast.pdb", dir, i), wwwcord, wwwmol, nil)

		xcord, xmol := BBHCap(c, m1, -1, v[len(v)-1], true)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dHOnlyC.pdb", dir, i), xcord, xmol, nil)
		xxcord, xxmol := BBHCap(c, m1, v[0], -1, false)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dHOnlyNPutMed.pdb", dir, i), xxcord, xxmol, nil)
		xxxcord, xxxmol := BBHCap(c, m1, v[0], -1, true)
		chem.PDBFileWrite(fmt.Sprintf("%s/results/cap%dHOnlyNPutLast.pdb", dir, i), xxxcord, xxxmol, nil)

	}

}
