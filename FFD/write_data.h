/* ----------------------------------------------------------------------------

	File:	        write.h

	Written by:   Wangda Zuo

	Task:         Declare subroutines in write.c              

---------------------------------------------------------------------------- */

int write_data(PARA_DATA *para, REAL **var, char *name);

void corners(PARA_DATA *para, REAL **var, REAL *psi);

int write_data1(PARA_DATA *para, REAL **var, char *name);
