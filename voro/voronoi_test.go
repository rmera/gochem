package voro

import (
	"fmt"
	"testing"

	chem "github.com/rmera/gochem"
	"github.com/rmera/gochem/solv"
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

func TTestFPlanes(t *testing.T) {
	mol, err := chem.PDBFileRead("../test/2c9vIOH.pdb", true)
	if err != nil {
		panic(err.Error())
	}
	res := []int{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 145, 146, 147, 148, 149, 150, 151, 152, 153}
	indexA := chem.Molecules2Atoms(mol, res, []string{"A"})
	indexB := chem.Molecules2Atoms(mol, res, []string{"F"})
	aindexes := make([]int, 0, len(indexA)+len(indexB))
	aindexes = append(aindexes, indexA...)
	aindexes = append(aindexes, indexB...)
	coord := mol.Coords[0]
	mol.FillVdw()
	fmt.Println("indexes", len(indexA), len(indexB))
	options := DefaultOptions()
	options.Subset = aindexes
	options.NoH = true
	cPlanes := ContactPlanes(coord, mol, options)
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
	fmt.Println("Atoms in Contact: ", len(ABConts), ABConts)
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
	res := []int{2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 145, 146, 147, 148, 149, 150, 151, 152, 153}
	indexA := chem.Molecules2Atoms(mol, res, []string{"A"})
	indexB := chem.Molecules2Atoms(mol, res, []string{"F"})
	aindexes := make([]int, 0, len(indexA)+len(indexB))
	aindexes = append(aindexes, indexA...)
	aindexes = append(aindexes, indexB...)
	coord := mol.Coords[0]
	mol.FillVdw()
	//	subcoord := v3.Zeros(len(aindexes))
	//	subcoord.SomeVecs(coord, aindexes) //all of them, in this case, but I'll keep this
	fmt.Println("indexes", len(indexA), len(indexB))
	options := DefaultOptions()
	options.Subset = aindexes
	options.NoH = false
	cPlanes := ContactPlanes(coord, mol, options)
	contacts := cPlanes.AllContacts()
	var ABConts [][2]int
	fmt.Println(len(contacts)) ///////////
	//	testatoms := indexA[2:6] ////////
	for _, v := range contacts {
		if (isInInt(indexA, v[0]) && isInInt(indexB, v[1])) || (isInInt(indexB, v[0]) && isInInt(indexA, v[1])) {
			ABConts = append(ABConts, v)
		}
	}
	var surf float64 = 0.00
	var vol float64
	var ctm []float64
	for _, v := range ABConts {

		//	if isInInt(indexA, v[0]) {
		ctm = PairContactAreaAndVolume(v[0], v[1], coord, mol, cPlanes)
		//	} else {
		//	ctm = PairContactAreaAndVolume(v[1], v[0], coord, cPlanes)

		//	}
		surf += ctm[0]
		vol += ctm[1]
		fmt.Println("Area and vol so far", surf, vol, "added now", ctm)
	}
	fmt.Println("Total contact area:", surf, "A^2. volume for the polyhedron associated to the 'A' chain part of the interface:", vol, "A^3")

}

/***************
**
**
**
**
**
**
**
**
**
****************/

func TestSolvAreas(t *testing.T) {
	fmt.Println("Will test solvated contact areas!!")
	mol, err := chem.PDBFileRead("../test/WTFull.pdb", true)
	if err != nil {
		panic(err.Error())
	}
	res := []int{151, 152, 153, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 145, 146, 147, 148, 149, 150, 151, 152, 153}
	indexA := chem.Molecules2Atoms(mol, res, []string{"A"})
	indexB := chem.Molecules2Atoms(mol, res, []string{"B"})
	aindexes := make([]int, 0, len(indexA)+len(indexB))
	aindexes = append(aindexes, indexA...)
	aindexes = append(aindexes, indexB...)
	coord := mol.Coords[0]
	mol.FillVdw()
	//solvation
	options := solv.DefaultOptions()
	options.Cpus(1)
	//options.COM(true)
	//options.End(20)
	solv := solv.DistRank(coord, mol, aindexes, []string{"SOL", "NA+", "CL-"}, options)
	var distanceCutoff float64 = 3.0
	solvdista := solv.UpToDistance(distanceCutoff)

	/////////////

	solvindexes := solvdista.AtomIDs(mol)
	//	fmt.Println(len(solvdista), solvindexes) //////////

	aindexes = append(aindexes, solvindexes...)
	redumol := chem.NewTopology(0, 1, nil)
	redumol.SomeAtoms(mol, aindexes)
	reduc := v3.Zeros(len(aindexes))
	reduc.SomeVecs(coord, aindexes)
	chem.PDBFileWrite("reduced.pdb", reduc, redumol, nil)
	//end solvation part
	//fmt.Println("indexes", len(indexA), len(indexB), len(solvindexes))
	scanoptions := DefaultOptions()
	scanoptions.Subset = aindexes
	cPlanes := ContactPlanes(coord, mol, scanoptions)
	contacts := cPlanes.AllContacts()
	var ABConts [][2]int
	var open []int
	var allConts []int
	//	fmt.Println(len(contacts)) ///////////
	//	testatoms := indexA[2:6] ////////
	for _, v := range contacts {
		if (isInInt(indexA, v[0]) && isInInt(indexB, v[1])) || (isInInt(indexB, v[0]) && isInInt(indexA, v[1])) {
			ABConts = append(ABConts, v)
		}
	}
	var surf float64 = 0.00
	var vol float64
	var ctm []float64
	for _, v := range ABConts {

		if isInInt(indexA, v[0]) {
			ctm = PairContactAreaAndVolume(v[0], v[1], coord, mol, cPlanes)
		} else {
			ctm = PairContactAreaAndVolume(v[1], v[0], coord, mol, cPlanes)

		}
		if ctm[0] == 0 {
			open = appendNotRepeatedInt(open, v[0])
			open = appendNotRepeatedInt(open, v[1])
		}
		surf += ctm[0]
		vol += ctm[1]
		fmt.Println("Area and vol so far", surf, vol, "added now", ctm, v)
	}
	if len(open) > 0 {
		openmol := chem.NewTopology(0, 1, nil)
		openmol.SomeAtoms(mol, open)
		openc := v3.Zeros(len(open))
		openc.SomeVecs(coord, open)
		chem.PDBFileWrite("open_polyhedra.pdb", openc, openmol, nil)
	}
	for _, v := range ABConts {
		allConts = appendNotRepeatedInt(allConts, v[0])
		allConts = appendNotRepeatedInt(allConts, v[1])
	}
	intermol := chem.NewTopology(0, 1, nil)
	intermol.SomeAtoms(mol, allConts)
	interc := v3.Zeros(len(allConts))
	interc.SomeVecs(coord, allConts)
	chem.PDBFileWrite("interface.pdb", interc, intermol, nil)

	fmt.Println("Total contact area:", surf, "A^2. volume for the polyhedron associated to the 'A' chain part of the interface:", vol, "A^3")

}