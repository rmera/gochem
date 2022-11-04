package voro

import (
	"fmt"
	"testing"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	//	"github.com/rmera/scu"
)

func TessssssRotation(t *testing.T) {
	mol, _ := chem.XYZFileRead("rota.xyz")
	coord := mol.Coords[0]
	axis := v3.Zeros(1)
	d := v3.Zeros(1)
	v1 := v3.Zeros(1)
	v2 := v3.Zeros(1)
	v2.Set(0, 1, 2.0) //any vector
	for _, angle := range []float64{10.0, 30.0, 90.0} {
		o := coord.VecView(1)
		d1 := coord.VecView(2)
		v1.Sub(d1, o)
		v2.Sub(v2, o)
		axis.Cross(v1, v2)
		ar := angle * chem.Deg2Rad
		d, _ = chem.RotateAbout(d1, o, axis, ar)
		coord.SetVecs(d, []int{2})
		name := fmt.Sprintf("rota_%3.1f.xyz", angle)
		chem.XYZFileWrite(name, coord, mol)
	}
}

func TTestVoronoi(t *testing.T) {
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

func TestFPlanes(t *testing.T) {
	mol, err := chem.PDBFileRead("../test/2c9vIOH.pdb", true)
	if err != nil {
		panic(err.Error())
	}
	res := []int{}
	for i := 1; i <= 153; i++ {
		res = append(res, i)
	}
	indexA := chem.Molecules2Atoms(mol, res, []string{"A"})
	indexB := chem.Molecules2Atoms(mol, res, []string{"F"})
	aindexes := make([]int, 0, len(indexA)+len(indexB))
	aindexes = append(aindexes, indexA...)
	aindexes = append(aindexes, indexB...)
	coord := mol.Coords[0]
	subcoord := v3.Zeros(len(aindexes))
	subcoord.SomeVecs(coord, aindexes) //all of them, in this case, but I'll keep this
	fmt.Println("indexes", len(indexA), len(indexB))
	cPlanes := ContactPlanes(subcoord, nil, nil)
	contacts := cPlanes.AllContacts()
	var ABConts []int
	fmt.Println(len(contacts)) ///////////
	//	testatoms := indexA[2:6] ////////
	for _, v := range contacts {
		if (isInInt(indexA, v[0]) && isInInt(indexB, v[1])) || (isInInt(indexB, v[0]) && isInInt(indexA, v[1])) {
			if !isInInt(ABConts, v[0]) {
				ABConts = append(ABConts, v[0])
			}
			if !isInInt(ABConts, v[1]) {
				ABConts = append(ABConts, v[1])
			}
		}
	}
	if len(ABConts) == 0 {
		t.Fatal("No contact found")
	}
	fmt.Println("Contacts: ", len(ABConts), ABConts)
	icoord := v3.Zeros(len(ABConts))
	icoord.SomeVecs(coord, ABConts)
	imol := chem.NewTopology(0, 1)
	imol.SomeAtoms(mol, ABConts)
	chem.PDBFileWrite("../test/inter.pdb", icoord, imol, nil)

}

func TTestAreas(t *testing.T) {
	mol, err := chem.PDBFileRead("../test/2c9vIOH.pdb", true)
	if err != nil {
		panic(err.Error())
	}
	const cutoff = 4
	mol.FillVdw()
	//res := []int{148, 149, 150, 151, 152, 153, 48, 49, 50, 51, 52}
	//	resB := []int{48, 49, 50, 51, 52}
	res := []int{}
	for i := 1; i <= 153; i++ {
		res = append(res, i)
	}
	indexA := chem.Molecules2Atoms(mol, res, []string{"A"})
	indexB := chem.Molecules2Atoms(mol, res, []string{"F"})
	coord := mol.Coords[0]
	planes := GetPlanes(coord, mol, cutoff, true)
	var ABConts [][]int
	//	testatoms := indexA[2:6] ////////
	for _, v := range indexA {
		for _, w := range indexB {
			angles := DefaultAngleScan() //this is the cutoff for inter-atom distances, so half of it is about right for atom-plane
			if planes.VdwContact(coord, mol, v, w, angles) {
				ABConts = append(ABConts, []int{v, w})
			}
		}
	}
	fmt.Println("Ready with A-B contacts")
	fmt.Println("Got", len(planes.ConfirmedContacts()), "A-B confirmed contacts")
	for i, v := range indexA {
		for _, w := range indexA[i:] {
			angles := DefaultAngleScan()                //this is the cutoff for inter-atom distances, so half of it is about right for atom-plane
			planes.VdwContact(coord, mol, v, w, angles) //we only want to "mark" the contact planes
		}
	}
	fmt.Println("Ready with A-A contacts")
	area := 0.0
	contactplanes := planes.ConfirmedContacts()
	fmt.Println("Got", len(contactplanes), "confirmed contacts")
	for _, v := range ABConts {
		area += PairContactArea(v[0], v[1], coord, contactplanes)
		fmt.Println("Area so far", area)
	}
	fmt.Println("Total contact area:", area, "A")

}
