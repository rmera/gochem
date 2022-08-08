#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
int get_coords(XDRFILE *fp, float *coordbuffer, float *box, int natoms); //Delcare like a good guy.
XDRFILE *openfile(char *name);
int read_natoms(char *name);
void xtc_close(XDRFILE *fp);
