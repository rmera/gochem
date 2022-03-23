package amberold

import (
	"flag"
	"fmt"

	//	"fmt"
	"os"
	"testing"

	chem "github.com/rmera/gochem"
	v3 "github.com/rmera/gochem/v3"
	//	"github.com/rmera/gochem/xtc"
	//	"github.com/rmera/scu"
	//	"gonum.org/v1/gonum/mat"
	//	"math"
	//	"sort"
	//	"strconv"
)

////use:  program [-skip=number -begin=number2] pdbfile trajname
func TestAmberold(Te *testing.T) {
	//The skip options
	skip := 0  // flag.Int("skip", 0, "How many frames to skip between reads.")
	begin := 1 //flag.Int("begin", 1, "The frame from where to start reading.")
	end := 100

	flag.Parse()
	//	println("SKIP", *skip, *begin, args) ///////////////////////////
	mol, err := chem.PDBFileRead("../../test/310K.pdb", false)
	if err != nil {
		Te.Error(err)

	}
	var traj chem.Traj
	fmt.Println("atoms:", mol.Len())
	trajname := "../../test/MDCuII_1_7000.crd"
	traj, err = New(trajname, mol.Len(), true) //false)
	if err != nil {
		Te.Error(err)
	}
	Coords := make([]*v3.Matrix, 0, 0)
	var coords *v3.Matrix
	lastread := -1
	for i := 0; i < end; i++ {
		if lastread < 0 || (i >= lastread+(skip) && i >= (begin)-1) {
			coords = v3.Zeros(traj.Len())
		}
		err := traj.Next(coords) //Obtain the next frame of the trajectory.
		if err != nil {
			_, ok := err.(chem.LastFrameError)
			if ok {
				break //We processed all frames and are ready, not a real error.

			} else {
				Te.Error(err)
			}
		}
		if (lastread >= 0 && i < lastread+(skip)) || i < (begin)-1 { //not so nice check for this twice
			continue
		}
		lastread = i
		Coords = append(Coords, coords)
		coords = nil // Not sure this works
	}
	pdbname := "../../test/MDCuII_1_7000.pdb"
	fout, err := os.Create(pdbname)
	if err != nil {
		Te.Error(err)
	}
	defer fout.Close()
	err = chem.MultiPDBWrite(fout, Coords, mol, nil)
	if err != nil {
		Te.Error(err)
	}
}
