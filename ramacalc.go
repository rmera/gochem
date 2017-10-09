package chem

import (
	"fmt"
	"math"
	"strings"
	"github.com/rmera/gochem/v3"
)


type RamaSet struct {
	Cprev   int
	N       int
	Ca      int
	C       int
	Npost   int
	MolID   int
	Molname string
}


// RamaCalc Obtains the values for the phi and psi dihedrals indicated in []Ramaset, for the
// structure M. The angles are in *degrees*.  It returns a slice of 2-element slices, one for the phi the next for the psi
// dihedral, a and an error or nil.
func RamaCalc(M *v3.Matrix, dihedrals []RamaSet) ([][]float64, error) {
	if M == nil || dihedrals == nil {
		return nil, CError{string(ErrNilData), []string{"RamaCalc"}}
	}
	r, _ := M.Dims()
	Rama := make([][]float64, 0, len(dihedrals))
	for _, j := range dihedrals {
		if j.Npost >= r {
			return nil, CError{"Data out of range", []string{"RamaCalc"}}
		}
		Cprev := M.VecView(j.Cprev)
		N := M.VecView(j.N)
		Ca := M.VecView(j.Ca)
		C := M.VecView(j.C)
		Npost := M.VecView(j.Npost)
		phi := Dihedral(Cprev, N, Ca, C)
		psi := Dihedral(N, Ca, C, Npost)
		temp := []float64{phi * (180 / math.Pi), psi * (180 / math.Pi)}
		Rama = append(Rama, temp)
	}
	return Rama, nil
}

// Filter the set of dihedral angles of a ramachandran plot by residue.(ex. only GLY, everything but GLY)
// The 3 letter code of the residues to be filtered in or out is in filterdata, whether they are filter in
// or out depends on shouldBePresent. It returns the filtered data and a slice containing the indexes in
// the new data of the residues in the old data, when they are included, or -1 when they are not included.
func RamaResidueFilter(dihedrals []RamaSet, filterdata []string, shouldBePresent bool) ([]RamaSet, []int) {
	RetList := make([]RamaSet, 0, 0)
	Index := make([]int, len(dihedrals))
	var added int
	for key, val := range dihedrals {
		isPresent := isInString(filterdata, val.Molname)
		if isPresent == shouldBePresent {
			RetList = append(RetList, val)
			Index[key] = added
			added++
		} else {
			Index[key] = -1
		}
	}
	return RetList, Index
}

// RamaList takes a molecule and returns a slice of RamaSet, which contains the
// indexes for each dihedral to be included in a Ramachandran plot. It gets the dihedral
// indices for all residues in the range resran. if resran has 2 elements defining the
// boundaries. Otherwise, returns dihedral lists for the residues included in
// resran. If resran has 2 elements and the last is -1, RamaList will
// get all the dihedral for residues from resran[0] to the end of the chain.
// It only obtain dihedral lists for residues belonging to a chain included in chains.
func RamaList(M Atomer, chains string, resran []int) ([]RamaSet, error) {
	RamaList := make([]RamaSet, 0, 0)
	if len(resran) == 2 {
		if resran[1] == -1 {
			resran[1] = 999999999 //should work!
		}
	}
	if M == nil {
		return nil, CError{"Nil data given", []string{"RamaList"}}
	}
	C := -1
	N := -1
	Ca := -1
	Cprev := -1
	Npost := -1
	chainprev := "NOTAVALIDCHAIN" //any non-valid chain name
	for num := 0; num < M.Len(); num++ {
		at := M.Atom(num)
		//First get the indexes we need. Change: If you give RamaList an empty string for "chains", it will
		//include all chains in the chem.Atomer.
		if strings.Contains(chains, string(at.Chain)) || at.Chain == " " || chains == "" {
			if at.Chain != chainprev {
				chainprev = at.Chain
				C = -1
				N = -1
				Ca = -1
				Cprev = -1
				Npost = -1
			}
			if at.Name == "C" && Cprev == -1 {
				Cprev = num
			}
			if at.Name == "N" && Cprev != -1 && N == -1 && at.MolID > M.Atom(Cprev).MolID {
				N = num
			}
			if at.Name == "C" && Cprev != -1 && at.MolID > M.Atom(Cprev).MolID {
				C = num
			}
			if at.Name == "CA" && Cprev != -1 && at.MolID > M.Atom(Cprev).MolID {
				Ca = num
			}
			if at.Name == "N" && Ca != -1 && at.MolID > M.Atom(Ca).MolID {
				Npost = num
			}
			//when we have them all, we save
			if Cprev != -1 && Ca != -1 && N != -1 && C != -1 && Npost != -1 {
				//We check that the residue ids are what they are supposed to be
				r1 := M.Atom(Cprev).MolID
				r2 := M.Atom(N).MolID
				r2a := M.Atom(Ca).MolID
				r2b := M.Atom(C).MolID
				r3 := M.Atom(Npost).MolID
				if (len(resran) == 2 && (r2 >= resran[0] && r2 <= resran[1])) || isInInt(resran, r2) {
					if r1 != r2-1 || r2 != r2a || r2a != r2b || r2b != r3-1 {
						return nil, CError{fmt.Sprintf("Incorrect backbone Cprev: %d N-1: %d CA: %d C: %d Npost-1: %d", r1, r2-1, r2a, r2b, r3-1), []string{"RamaList"},}
					}
					temp := RamaSet{Cprev, N, Ca, C, Npost, r2, M.Atom(Ca).Molname}
					RamaList = append(RamaList, temp)
				}
				N = Npost
				Ca = -1
				Cprev = C
				C = -1
				Npost = -1
			}
		}
	}
	//	fmt.Println("Rama",Rama, "failed", failed)
	return RamaList, nil
}
