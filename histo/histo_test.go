package histo

import (
	"encoding/json"
	"fmt"
	"testing"
)

func TestHistoIO(Te *testing.T) {
	fmt.Println("Histogram JSON output test!")
	M := NewMatrix(3, 3, []float64{0, 1, 2, 3, 4, 8})
	M.Fill()
	rawdata := []float64{1, 6, 3, 2, 4, 5, 7, 6, 3.5, 3, 5, 1, 1, 0, 0, 5, 8, 1, 2, 3, 44, 3, 7, 3, 1, 3, 5, 32, 1}
	M.NewHisto(0, 1, nil, rawdata)
	v := M.View(0, 1)
	fmt.Println(v.String())
	j, err := json.Marshal(M)
	fmt.Println("JSON:", string(j), err)
	M2 := new(Matrix)
	json.Unmarshal(j, M2)
	fmt.Printf("%v\n", M2)

}
