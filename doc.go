/*
 * doc.go, part of gochem.
 *
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * Gochem is developed at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland.
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

/*Package chem is the main package of the goChem library. It provides atom and molecule structures, facilities for reading and writing some
files used in computational chemistry and some functions for geometric manipulations and shape
indicators.



	**goChem Capabilities**


    Reads/writes PDB and XYZ files.

    Reads XTC and DCD trajectory files, both sequentially and concurrently.

    Superimposes molecules (especially adequate for non-proteins since doesn't
	use sequence information). The user specify what atoms to use for the
	superimposing transformation calculation. Then all the atoms will be
	superimposed accordingly. Thus, non-identical molecules can be superimposed.

    Calculates RMSD between sets of coordinates.

    Allows to select atoms and coordinates by using a go slice of indexes.

    Allows to replace selected coordinates for a new set.

    Calculates moment tensor and elipsoid of inertia--related properties

    The Molecule object implements the sort.Interface interface, so atoms
	can easily be sorted by b-factors.

    Axis manipulation.
        Align a vector with the Z axis.
        Rotate around the Z axis until the xy projection of a vector becomes
		the Y axis.
        Rotate a sub-group of atoms in a molecule using any 2 coordinates as
		the rotation axis.

    The latter is implemented using Clifford algebra and, as a legacy version, Euler
	angles and rotation matrices (math for the Cliffor algebra implementation by Doc.
	Dr. Janne Pesonen). The Clifford algebra implementation is concurrent. In general,
	Clifford algebra is mathematically better behaving than Euler angles, which are
	not defined for certain rotations.

    Calculates and draws Ramachandran plots (uses the Plotinum library). for an
	aminoacidic chain or a subset of it.

    Generates input for, run and recover results from QM calculations with Turbomole,
	Orca and MOPAC (which must be obtained independently from their respective
	distributors). Interfacing gochem to other QM codes is fairly simple.

    goChem data can be JSON encoded and transfered in such a way that PyMOL
	(http://www.pymol.org plugins can send and received data from/to goChem programs)
	thus enabling the build of PyMOL plugins using goChem.



goChem implements its own matrix type for coordinates, VecMatrix,based in github.com/gonum/matrix.

Currently, each row of a VecMatrix represents one point in space. As this could change if
github.com/gonum/matrix changes, we recomend prefering the Vec* methods over the Row* methods
when manipulating a VecMatrix. Vec* methods will change from row to column following the upstream
Gonum library.*/
package chem
