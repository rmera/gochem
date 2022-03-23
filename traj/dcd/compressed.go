/*
 * dcd.go, part of gochem
 *
 * Copyright 2012 Raul Mera Adasme <rmera_changeforat_chem-dot-helsinki-dot-fi>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License  as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

package dcd

import (
	"bufio"
	"compress/flate"
	"compress/lzw"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

const (
	lzwOrder        = lzw.MSB
	lzwLitwidth int = 8
)

//prepSource takes a filename and format string, opens the file and returns an object that will
// read data from the file, either 'as is' or decompressing first, it depending on the format string.
//If the format string is empty, it will try to deduce it form the file extension. File extensions supported are
//.dcd (non-compressed dcd), .gz (deflate) and .lzw. If the format string is empty and the extension doesn't
//match any supported type, a message will be logged and the non-compressed dcd format will be assumed.
//thus, prepSource only returns an error if the file can't be opened.
func (D *DCDObj) prepSource(fname string, format string) (io.ReadCloser, error) {
	var err error
	var fk string
	var ret io.ReadCloser
	if format == "" {
		temp := strings.Split(fname, ".")
		fk = strings.ToLower(temp[len(temp)-1])
	} else {
		fk = format
	}
	D.filename = fname
	D.fhandle, err = os.Open(fname)
	if err != nil {
		return nil, Error{err.Error(), D.filename, []string{"os.Open", "prepSource"}, true}
	}
	reader := bufio.NewReader(D.fhandle)

	switch fk {

	case "dcd":
		ret = D.fhandle

	case "lzw":
		ret = lzw.NewReader(reader, lzwOrder, lzwLitwidth)
	case "gz":
		ret = flate.NewReader(reader)
	default:
		//if it's not a plain DCD, you'll get an error later.
		log.Printf("Format string %s not supported. %s will be assumed to be a plain DCD file", fk, D.filename)
		D.dcd = D.fhandle
	}
	return ret, nil

}

//parselevel parses a string into a compression extension and
//a compression level.
//The string needs to be  in the format: gz_level (ejemplo: gz_-2)
//if any part of the procedure fails, we will assume that there is no
//"level" and we return the whole file name, and a default value for
//level (-2, corresponding to the poorest, but faster, compression.
//Note that the compression level obtained only makes sense in the
//context of deflate (zlib or gzip) compression.
func parseLevel(f string) (string, int) {
	data := strings.Split(f, "_")
	if len(data) != 2 || data[0] != strings.ToLower("gz") {
		return f, -2
	}
	level, err := strconv.Atoi(data[1])
	if err != nil {
		return f, -2
	}
	if level < -2 || level > 9 {
		if data[0] == "gz" {
			return data[0], -2
		} else {
			return f, -2
		}
	}
	return data[0], level
}

//As you can see I wrote the whole function just to be reminded that it can not work, because
//DCD for _some_ reason, requires the number of frames at the begining, and so you need to
//move back and forth through the file every time you write a new frame. ^_^ We really
//need less crazy formats for MD.

//prepTarget takes fname, the name of a DCD trajectory file and format, a string
//that is either empty or indicates a compression format.
//It creates the file and returns a io.WriteCloser that will write data, crude or
//compressed, with deflate (gzip/zlib) or with lzw, depending on the file extension.
//(.gz for deflate, .lzw for lzw, any other extension for a non compressed dcd)
//The file extension may contain information about the compression level,
//which applies only to deflate. If so, the requested compression level will be used.
/*
prepTarget(D *DCDWObj, fname string) (io.WriteCloser, error) {
	var err error
	var fk string
	var tmpfk string
	var ret io.WriteCloser
	temp := strings.Split(fname, ".")
	tmpfk = strings.ToLower(temp[len(temp)-1])
	fk, level := parseLevel(tmpfk)
	D.filename = fname
	D.fhandle, err = os.Create(fname)
	if err != nil {
		return nil, Error{err.Error(), D.filename, []string{"os.Open", "prepTarget"}, true}
	}
	writer := bufio.NewWriter(D.fhandle)

	switch fk {

	case "dcd":
		ret = D.fhandle

	case "lzw":
		ret = lzw.NewWriter(writer, lzwOrder, lzwLitwidth)
	case "gz":
		ret, err = flate.NewWriter(writer, level)
	default:
		log.Printf("Format string %s not supported. %s will be written to be a plain DCD file", fk, D.filename)
		D.dcd = D.fhandle
	}
	return ret, nil

}
**/
