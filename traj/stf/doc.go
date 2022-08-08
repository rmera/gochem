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
//can be easily implemented in other programing languages / for other libraries or programs.

/******************** Format Specification   ***************************************************

The substrings "*" and "**" can't appear anywhere in the file except in the parts indicated.

An STF file has the extension stf, and it is compressed with z-standard (zstd)

A STF file may only contain ASCII symbols.

A STF file has a "header" starting in the first line, and ending with a line that contains "**" followed by one or more spaces, and the number of atoms per frame.
Each line of the header must be a pair key=value.

After that, the file has one line per atom, per frame. Each line 3 numbers (x y z coordinates, in A). The precision is not specified.
each frame ends with a line starting with the character (no whitespaces before) "*", optionally followed by a whitespace and 9 floating-point numbers separated by spaces (precision unspecified). If present, these number correspond to the vectors defining the simulation box, in A.

The "**" sequence may not be used anywhere in the file, except in the line described above.

***************************************************************************************************/

package stf
