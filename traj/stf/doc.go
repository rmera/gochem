/*
 * doc.go, part of gochem.
 *
 * Copyright 2021 Raul Mera <rauldotmeraatusachdotcl>
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
 *
 * /

//package stf implements the simple trajectory format, an internal trajectory format for goChem.
//stf aimst to produce reasonably small files and to be very easy to read and write, so readers/writers
//can be easily implemented in other programing languages / for other libraries or programs, while
//also being reasonably fast to write and, especially, to read.

//By my test, stf files are ~60% larger than the equivalent XTC file at equivalent precision (and about
30% larger at the lower default stf precision), while being about half the size of a an equivalent DCD file


/******************** Format Specification   ***************************************************


An STF file has the extension stf, and it is compressed with z-standard (zstd).

A STF file may only contain ASCII symbols.

A STF file has a "header" starting in the first line, and ending with a line that starts with the
characters "**" followed by one or more spaces, and the number of atoms per frame.

Each line of the header must be a pair key=value. It a trajectory includes topology, this may be
included in the header with the key "topology" and a jsons string, describing the topology as an
array of atoms, as defined in github.com/rmera/gochem/chemjson, as a value. Implementations are
not required to support reading/writing of topologies. The precision (an integer greater than 0,
see below) must be included in the header, with the corresponding key "prec". For example,
a 'precision' line could be:

prec=2

Note that the header must have at least one line, marking the precison. The 'default' precision
to be used if the user does not supply one is left to the implementation, though we do note that
the implementation in this package uses a default precision of 1.

After the header, the file has one line per atom, per frame. Each line contains  3 numbers,
corresponding to the x y and z cartesian coordinates, respectively, and nothing more. Each
of these 3 number contains the respective coordinate in Angstrom, multiplied by 10 to the
power of (precision), where precision is 1, or whatever positive integer given in the header
(see above) and rounded to make it an integer. Since the default precision is
1, the numbers go multiplied by 10, unless a different value is given in the header, and rounded
to an integer.

Each frame ends with a line starting with the character "*" (no whitespaces before) , optionally
followed by: one or more whitespace and 9 floating-point numbers separated by spaces (precision
unspecified). If present, these number correspond to the vectors defining the simulation box, in Angstrom

The "**" sequence may only be used as a header termination, as described above and can not appear
anywhere else in the file.

Z-standard allows several 'levels' of compression that represent different tradeoffs between
compression/decompression speed and final file size. The naming of those levels is implementation
dependent. The level used in stf is left to the implementation, but at least a default level must
be provided. If an implementation does _not_ give the user an option to control the level used
(we recommend, but do not require, that such option is given), we request lack of option to be
clearly documented.

***************************************************************************************************/

package stf
