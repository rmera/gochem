package align

import "runtime"

//Options contains various options for the LOVOnMostRigid and RMSDTraj functions
type Options struct {
	begin int
	skip  int
	cpus  int
	//The following are ignored by the RMSDTraj function
	atomNames []string
	chains    []string
	//	nMostRigid   int
	writeTraj    string  //the name of a file to where the aligned trajectory will be written. Nothing will be written if empty.
	lessThanRMSD float64 //instead of using a N for the most rigid, selects all residues with RMSD (the square root of the RMSD) < LessThanRMSD, in A. If >0, this overrrides NMostRigid
	minimumN     int     //The smallest N we are willing to have. Only valid if LessThanRMSD is in use.
}

//DefaultOptions return reasonable options for atomistic trajectories.
//It prepares a superposition of alpha carbons (CA) with all logical CPUs,
//trying to use for the superposition all CAs with RMSD lower than 1.0 A.
func DefaultOptions() *Options {
	r := new(Options)
	r.cpus = runtime.NumCPU()
	//all available CPUs
	r.atomNames = []string{"CA"}
	r.chains = nil
	//	r.nMostRigid = -1
	r.lessThanRMSD = 1.0
	r.minimumN = 10 //just a reasonable value.
	return r
}

//DefaultCGOptions returns reasonable options for Martini trajectories.
func DefaultCGOptions() *Options {
	r := DefaultOptions()
	r.atomNames = []string{"BB"}
	//	r.Skip = 1000 //seems reasonable for CG, those trajectories are _long_.
	return r

}

//Sets O.N to represent the perc percent of the residues
//in the sequence. It also requires seqlen, the total number of
//residues in the system
//func (O *Options) SetRigidPercent(perc int, seqlen int) {
//	frac := float64(perc) / 100
//	O.nMostRigid = int(frac * float64(seqlen))
//
//}

//Returns the value of the
//func (O *Options) NMostRigid(N ...int) int {
//	if len(N) < 0 || N[0] <= 0 {
//		return O.nMostRigid
//	}
//	O.nMostRigid = N[0]
//	return N[0]
//
//}

//Returns the value of the first frame of the sequence to use,
//and sets it to a new value, if given.
func (O *Options) Begin(n ...int) int {
	if len(n) > 0 || n[0] >= 0 {
		O.begin = n[0]
	}
	return O.begin
}

//Returns the value of the skipped frames between reads,
//and sets it to a new value, if given.
func (O *Options) Skip(n ...int) int {
	if len(n) > 0 || n[0] >= 0 {
		O.skip = n[0]
	}
	return O.skip
}

//Returns the number of gorutines to be used,
//and sets it to a new value, if given.
func (O *Options) Cpus(n ...int) int {
	if len(n) > 0 || n[0] > 0 {
		O.cpus = n[0]
	}
	return O.cpus
}

//
// The following options are ignored by the RMSDTraj function
//

//Returns the minimum RMSD for an atom to be consiered
//as part of the alignment, and sets it to a new value,
//if given.
//if this flag is set to 0 or less, the regular LOVO
//convergency criterion is used.
func (O *Options) LessThanRMSD(rmsd ...float64) float64 {
	if len(rmsd) > 0 {
		O.lessThanRMSD = rmsd[0]
	}
	return O.lessThanRMSD

}

//Returns the smallest acceptable number of atoms
//to be returned by the alinment procedure. Only used
//if the LessThanRMSD is active.
//and sets it to a new value, if given.
func (O *Options) MinimumN(n ...int) int {
	if len(n) > 0 && n[0] > 0 {
		O.minimumN = n[0]
	}
	return O.minimumN
}

//Returns the atom names consiered or the alignment
//and sets them to new values, if those are given.
func (O *Options) AtomNames(names ...[]string) []string {
	if len(names) > 0 && len(names[0]) > 0 {
		O.atomNames = names[0]
	}
	return O.atomNames
}

//Returns the atom names consiered or the alignment
//and sets them to new values, if those are given.
func (O *Options) Chains(chains ...[]string) []string {
	if len(chains) > 0 {
		O.chains = chains[0]
	}
	return O.chains
}

//Returns the name of the files where the aligned trajectory will be written, and sets it to a
//new value, if given. No trajectory file will be given if this value is set to an empty string.
func (O *Options) TrajName(name ...string) string {
	if len(name) > 0 && name[0] != "" {
		O.writeTraj = name[0]
	}
	return O.writeTraj
}
