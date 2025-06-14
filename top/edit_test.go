package top

import (
	"fmt"
	"testing"
)

func TestEdit(Te *testing.T) {
	m := AddOrDelAtomFunctions(true, 6, 1)
	res, err := FuncApplier("tops/Protein_Mb.itp", m)
	if err != nil {
		Te.Error(err)
	}
	fmt.Print(res)
}
