package prot

import (
	"fmt"
	"os"
	"slices"
	"strings"

	chem "github.com/rmera/gochem"
)

var maxProtons = map[string][]string{
	"LYS": []string{"HZ1", "HZ2", "HZ3"},
	"ARG": []string{"HE", "HH11", "HH12", "HH21", "HH22"},
	"ASP": []string{"HD2"},
	"GLU": []string{"HE2"},
}

// Checks if the residue is in it's common protonation state.
// first is the index of the first atom of the residue to be checked. Or any
// atom of the residue that comes before the labile protons that we'll be looking for.
func checkChargeRes(mol chem.Atomer, first int) int {
	pos := 0

	f := mol.Atom(first)
	MID := f.MolID
	name := f.MolName
	check := maxProtons[name]
	for i := first; i < mol.Len(); i++ {
		at := mol.Atom(i)
		if at.MolID != MID {
			break
		}
		if slices.Contains(check, at.Name) {
			pos++
		}
	}
	if pos == len(check) && slices.Contains([]string{"LYS", "ARG"}, name) {
		return 0
	} else if pos == len(check) {
		return +1 //means we have a protonated acid residue
	}
	return -1 //we have a deprotonated basic residue

}

// Returns the charge of the protein. In principle, it assumes all basic residues (but histidine) are protonated
// and all acid residues are deprotonated, but it can check that the correct hydrogens are present/absent acoording to AMBER
// atom names, if the carefulCheck flag is given and true.
func Charge(mol chem.Atomer, atomname string, carefulCheck ...bool) int {

	var chargedres []int //IDs of thbe alpha carbon of every potentially charge dresidue except for histidine.
	charge := 0
	readinghis := false
	HisHs := 0
	IDs := make([]int, 0, 2)
	for i := 0; i < mol.Len(); i++ {
		curr := mol.Atom(i)
		if !slices.Contains(IDs, curr.MolID) {
			IDs = append(IDs, curr.MolID)
		}
		if curr.Name == atomname && (curr.MolName == "GLU" || curr.MolName == "ASP") {
			chargedres = append(chargedres, i)
			charge -= 1
		}
		//I added HIP in case the user is nice enough to point out the charged histidines.
		//We try to detect them otherwise (see below)
		if curr.Name == atomname && (curr.MolName == "ARG" || curr.MolName == "LYS" || curr.MolName == "HIP") {
			if curr.MolName != "HIP" {
				chargedres = append(chargedres, i)
			}

			charge += 1
		}
		if curr.MolName == "HIS" {
			readinghis = true
		}
		//The last condition i == mol.Len()-1 is new. Before that, if you happened to
		//have a charges histidine as last residue, you didn't get its charge.
		if readinghis && curr.Name == "HD1" {
			HisHs += 1
		}
		if readinghis && curr.Name == "HE2" {
			HisHs += 1
		}
		if readinghis && HisHs == 2 { //We have a doubly-protonated histidine
			HisHs = 0 //reset the counter
			charge += 1
		}

	}
	if len(carefulCheck) > 0 && carefulCheck[0] {
		for _, v := range chargedres {
			charge += checkChargeRes(mol, v)
		}
	}

	return charge
}

// Returns the heavy atom bonded to the given hydrogen in one of the standard
// 20 aminoacids (I don't have selenocysteine yet)
func HeavyAtomBoundToH(aaname, atname string) (string, error) {
	if len(aaname) < 3 {
		return "", fmt.Errorf("")
	}
	//Unfortunately, I need to assume you have the environment variable set.
	fname := os.Getenv("GOCHEM") + "/prot/oaa/" + strings.ToUpper(aaname) + ".pdb"
	aa, err := chem.PDBFileRead(fname, true)
	if err != nil {
		return "", fmt.Errorf("Couldn't open file for residue %s: %w", aaname, err)
	}
	err = aa.AssignBonds(aa.Coords[0])
	if err != nil {
		return "", fmt.Errorf("Couldn't assign bonds to residue %s: %w", aaname, err)
	}

	var at *chem.Atom
	for i := 0; i < aa.Len(); i++ {
		at = aa.Atom(i)
		if at.Name == atname {
			break
		}
		at = nil
	}
	if at == nil {
		return "", fmt.Errorf("Couldn't find atom %s", atname)

	}
	if len(at.Bonds) > 1 {
		return "", fmt.Errorf(" Atom %s Doesn't seem to be an hydrogen. %d bonds were detected for it", aaname, len(at.Bonds))
	}
	at2 := at.Bonds[0].Cross(at)
	return at2.Name, nil
}
