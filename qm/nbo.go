package qm

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"

	chem "github.com/rmera/gochem"
)

// just for brevity
var con func(string, string) bool = strings.Contains
var Err func(string, ...any) error = fmt.Errorf
var trim func(string) string = strings.TrimSpace

type NBOData struct {
	Charges []*Charge
	EConfs  []*EConf
	Steric  []*Steric
	Delocs  []*Deloc
}

// Obtains NBO properties from the output file. Properties to obtain are:
// Natural Charges, delocalizations and steric energies.
func NBOProps(inpf io.Reader) (*NBOData, error) {
	inp := bufio.NewReader(inpf)
	ret := new(NBOData)
	var nbo bool
	var err error
	var l string
	for l, err = inp.ReadString('\n'); err == nil; l, err = inp.ReadString('\n') {
		if con(l, "N A T U R A L   B O N D   O R B I T A L   A N A L Y S I S") {
			nbo = true
			continue
		}
		if !nbo {
			continue
		}
		if con(l, "Total disjoint NLMO steric exchange energy from pairwise sum") {
			break
		}
		if con(l, "Summary of Natural Population Analysis:") {
			ret, inp, err = SectionReader(inp, npa, ret)
			if err != nil {
				return nil, Err("NPA Error %w", err)
			}
		}
		if con(l, "Atom No         Natural Electron Configuration") {
			ret, inp, err = SectionReader(inp, econf, ret)
			if err != nil {
				return nil, Err("Natural Configuration Error %w", err)
			}
		}
		if con(l, "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS") {
			ret, inp, err = SectionReader(inp, delocs, ret)
			if err != nil {
				return nil, Err("Delocalizations error %w", err)
			}

		}
		if con(l, "NBO/NLMO STERIC ANALYSIS:") {
			ret, inp, err = SectionReader(inp, steric, ret)
			if err != nil {
				return nil, Err("Steric error %w", err)
			}

		}
	}
	if errors.Is(err, io.EOF) {
		err = nil
	}
	return ret, nil
}

type section struct {
	ini    func(string) bool
	end    func(string) bool
	reader func(string, *NBOData) error
}

func print1or2(s []int) string {
	r := make([]string, 0, len(s))
	for _, v := range s {
		r = append(r, fmt.Sprintf("%d", v))
	}
	return strings.Join(r, "-")
}

type Deloc struct {
	Donors    []int //could be 1 or 2
	Acceptors []int //same
	DonorType string
	AccType   string
	Fock      float64
	E2        float64
	DeltaE    float64
}

func (D *Deloc) String() string {
	return fmt.Sprintf("D: %s %s | A: %s %s | E2: %4.1f F: %4.1f DE %4.1f", D.DonorType, print1or2(D.Donors), D.AccType, print1or2(D.Acceptors), D.E2, D.Fock, D.DeltaE)
}

func orbids(bondtype, line string, i, j, k, l int) ([]int, error) {
	don, err := strconv.Atoi(trim(line[i:j]))
	if err != nil {
		return nil, err
	}
	r := []int{don}
	if con(bondtype, "BD") {
		don, err := strconv.Atoi(trim(line[k:l]))
		if err != nil {
			return nil, err
		}
		r = append(r, don)
	}
	return r, nil

}

var firstorbsindexes [4]int = [4]int{16, 19, 22, 25}
var secorbsindexes [4]int = [4]int{44, 47, 51, 54}

func delocsreader(l string, nbo *NBOData) error {
	if l == "" || l == "\n" || con(l, "unit") {
		return nil
	}
	d := new(Deloc)
	d.DonorType = l[7:10]
	if d.DonorType[2] == ' ' {
		d.DonorType = d.DonorType[:2]
	}
	d.AccType = l[35:38]
	if d.AccType[2] == ' ' {
		d.AccType = d.AccType[:2]
	}
	var err error
	fi := firstorbsindexes[:]
	d.Donors, err = orbids(d.DonorType, l, fi[0], fi[1], fi[2], fi[3]) //the numbers are the positions of the data we want in the line.
	if err != nil {
		return err
	}
	s := secorbsindexes[:]
	d.Acceptors, err = orbids(d.AccType, l, s[0], s[1], s[2], s[3])
	if err != nil {
		return err
	}
	d.E2, err = strconv.ParseFloat(trim(l[59:63]), 64)
	if err != nil {
		return err
	}
	d.DeltaE, err = strconv.ParseFloat(trim(l[67:71]), 64)
	if err != nil {
		return err
	}
	d.DeltaE *= chem.H2Kcal
	d.Fock, err = strconv.ParseFloat(trim(l[74:]), 64)
	if err != nil {
		return err
	}
	d.Fock *= chem.H2Kcal
	nbo.Delocs = append(nbo.Delocs, d)

	return nil
}

var delocs section = section{ini: func(l string) bool {
	return con(l, "===============================================================================")
},
	end: func(l string) bool {
		return con(l, "NATURAL BOND ORBITALS (Summary):")
	},
	reader: delocsreader,
}

type Steric struct {
	Type   string
	OrbIDs []int
	E      float64
}

func (S *Steric) String() string {
	return fmt.Sprintf("Orb: %s %s E: %4.1f", S.Type, print1or2(S.OrbIDs), S.E)
}

var steric section = section{ini: func(l string) bool {
	return con(l, "Occupied NLMO contributions dE(i) (kcal/mol) to total steric exchange energy")
},
	end: func(l string) bool { return con(l, "-------------------------------------------------") },
	reader: func(l string, n *NBOData) error {
		if l == "" || con(l, "unit") || l == "\n" {
			return nil
		}
		var err error
		s := new(Steric)
		s.Type = trim(l[7:9])
		fi := firstorbsindexes[:]
		s.OrbIDs, err = orbids(s.Type, l, fi[0], fi[1], fi[2], fi[3])
		if err != nil {
			return err
		}
		s.E, err = strconv.ParseFloat(trim(l[35:41]), 64)
		if err != nil {
			return err
		}
		n.Steric = append(n.Steric, s)
		return nil
	},
}

type Charge struct {
	Symbol string
	ID     int
	Charge float64
}

func (C *Charge) String() string {
	return fmt.Sprintf("%s %d Charge: %3.1f", C.Symbol, C.ID, C.Charge)
}

var npa section = section{ini: func(l string) bool {
	return con(l, "--------------------------------------------------------------------")
},
	end: func(l string) bool {
		return con(l, "====================================================================")
	},
	reader: func(l string, n *NBOData) error {
		var err error
		if n == nil {
			return fmt.Errorf("Given nil NBOData!")
		}
		if n.Charges == nil {
			n.Charges = make([]*Charge, 0, 10)
		}
		f := strings.Fields(l)
		ch := new(Charge)
		ch.Symbol = f[0]
		ch.ID, err = strconv.Atoi(f[1])
		if err != nil {
			return err
		}
		ch.Charge, err = strconv.ParseFloat(f[2], 64)
		if err != nil {
			return err
		}
		n.Charges = append(n.Charges, ch)
		return err
	},
}

type EConf struct {
	Symbol string
	ID     int
	Conf   string
}

func (C *EConf) String() string {
	return fmt.Sprintf("%s %d %s", C.Symbol, C.ID, C.Conf)
}

var econf section = section{ini: func(l string) bool {
	return con(l, "----------------------------------------------------------------------------")
},
	end: func(l string) bool {
		if l == "" || l == "\n" {
			return true
		}
		return false
	},
	reader: func(l string, n *NBOData) error {
		var err error
		if n == nil {
			return fmt.Errorf("Given nil NBOData!")
		}
		if n.EConfs == nil {
			n.EConfs = make([]*EConf, 0, 10)
		}
		f := strings.Fields(l)
		c := new(EConf)
		c.Symbol = f[0]
		c.ID, err = strconv.Atoi(f[1])
		if err != nil {
			return err
		}
		c.Conf = trim(strings.Join(f[2:], ""))
		n.EConfs = append(n.EConfs, c)
		return err
	},
}

func SectionReader(inp *bufio.Reader, sec section, nbo *NBOData) (*NBOData, *bufio.Reader, error) {
	var reading bool
	var err error
	for l, err := inp.ReadString('\n'); err == nil; l, err = inp.ReadString('\n') {
		if sec.ini(l) {
			reading = true
			continue
		}
		if !reading {
			continue
		}
		if sec.end(l) {
			break
		}
		err := sec.reader(l, nbo)
		if err != nil {
			return nbo, nil, err
		}
	}
	return nbo, inp, err

}
