/*
 * xtc.c
 * 
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License  as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h"
#include "xdrfile_xtc.h"

int get_coords(XDRFILE *fp, float *coordbuffer, float *box, int natoms); 
XDRFILE *openfile(char *name);
int read_natoms(char *name);
void xtc_close(XDRFILE *fp);


/*
int main(int argc, char **argv){
	int first=1;
	int skip=2;
	int end=6;
	int i=0;  //iterator
	char name[20]="test.xtc";
	XDRFILE *fp;
	int natoms=0;
	float *coords;
	read_xtc_natoms(name,&natoms);
	if (natoms==0){
		return 1;
		}
	printf("read!, %d\n",natoms);
	fp=xdrfile_open(name,"r");
	if (fp==NULL){
		return 1;
		}
	for(i=first;i<=end;i++){
		coords=get_coords(fp,natoms);                
		}
	
	}
*/


XDRFILE *openfile(char *name){
	XDRFILE *fp;
	fp=xdrfile_open(name,"r");
	return fp; /*NULL if something goes wrong.*/
	}

int read_natoms(char *name){
	int natoms=0;
	read_xtc_natoms(name,&natoms);
	return natoms; /*should return 0 if something goes wrong.*/
	}
	
void xtc_close(XDRFILE *fp){
	xdrfile_close(fp);
	}

/*get_coords is meant to be called from Go functions.
It takes an XDRFILE pointer (fp), which can be treated just as an opaque 
* object in Go, the number of atoms to be read (natoms)
* and a *float pointer to a buffer with enough space to 
* save all the coordinates in a frame. This buffer should be allocated
* from Go, so it is garbage collected.
* The function returns 0 on success and other number depending on the 
* error encountered*/
int get_coords(XDRFILE *fp, float *coordbuffer,float *box, int natoms){
	int step=0;
	int success=0;
	rvec *coords;
	rvec *grobox;
	//float box[3][3];
	float prec=0;
	float time=0;
	coords=(rvec *)coordbuffer; /*This is the type the GROMACS function wants*/
	grobox=(rvec *)box;
	if (coords==NULL){
		return 1;
		}
	success=read_xtc(fp,natoms,&step,&time,grobox,
                        coords,&prec);
	return success;
	}

