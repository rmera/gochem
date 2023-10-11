package voro

import (
	"fmt"
	"testing"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	//	"github.com/rmera/scu"
)

func TestVoronoi(t *testing.T) {
	mol, err := chem.PDBFileRead("../test/2c9vIOH.pdb", true)
	if err != nil {
		panic(err.Error())
	}
	const cutoff = 4
	mol.FillVdw()
	//	resA := []int{148, 149, 150, 151, 152, 153}
	//	resB := []int{48, 49, 50, 51, 52}
	res := []int{}
	for i := 1; i <= 153; i++ {
		res = append(res, i)
	}
	indexA := chem.Molecules2Atoms(mol, res, []string{"A"})
	indexB := chem.Molecules2Atoms(mol, res, []string{"F"})
	coord := mol.Coords[0]
	planes := GetPlanes(coord, mol, cutoff, true)
	var ABConts []int
	//	testatoms := indexA[2:6] ////////
	for _, v := range indexA {
		for _, w := range indexB {
			angles := DefaultAngleScan() //this is the cutoff for inter-atom distances, so half of it is about right for atom-plane
			if planes.VdwContact(coord, mol, v, w, angles) {
				ABConts = append(ABConts, v, w)
			}
		}
	}
	if len(ABConts) == 0 {
		t.Fatal("No contact found")
	}
	fmt.Println("Contacts: ", len(ABConts))
	icoord := v3.Zeros(len(ABConts))
	icoord.SomeVecs(coord, ABConts)
	imol := chem.NewTopology(0, 1)
	imol.SomeAtoms(mol, ABConts)
	chem.PDBFileWrite("../test/inter.pdb", icoord, imol, nil)
}
