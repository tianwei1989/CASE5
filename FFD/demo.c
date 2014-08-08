/*======================================================================
 demo.c
=======================================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <glut.h>
#include <time.h>

#include "data_structure.h"
#include "solver.h"
#include "write_data.h"
#include "init.h"
#include "boundary.h"

#include "grid.h"
#include "timing.h"
#include "read_data.h"


/* global variables */
static REAL dt, diff, visc;
static REAL force, source;
static int screen;

REAL **var;
REAL *x, *y, *z, *gx, *gy, *gz;
REAL *u, *v, *w, *u_s, *v_s, *w_s, *u_mean, *v_mean, *w_mean;
REAL *dens, *dens_s, *temp, *temp_s, *temp_mean, *p, *my_div, *pp;
REAL *tmp1, *tmp2, *tmp3;
REAL *ap, *an, *as, *aw, *ae, *b, *ab, *af, *ap0;
REAL *flagp, *flagu, *flagv, *flagw;

static GEOM_DATA geom;
static PROB_DATA prob;
static TIME_DATA mytime;
static OUTP_DATA outp1;
static BC_DATA bc;
static SOLV_DATA solv;
static PARA_DATA para;


clock_t start, end;

/******************************************************************************
| allocate data
******************************************************************************/
int allocate_data ( void )
{
	int size = (geom.imax+2) * (geom.jmax+2) * (geom.kmax+2);

	x			    = (REAL *) malloc ( size*sizeof(REAL) );
  y			    = (REAL *) malloc ( size*sizeof(REAL) );
  z			    = (REAL *) malloc ( size*sizeof(REAL) );
  u			    = (REAL *) malloc ( size*sizeof(REAL) );
	v			    = (REAL *) malloc ( size*sizeof(REAL) );
  w			    = (REAL *) malloc ( size*sizeof(REAL) );
	u_s		    = (REAL *) malloc ( size*sizeof(REAL) );
	v_s		    = (REAL *) malloc ( size*sizeof(REAL) );
  w_s		    = (REAL *) malloc ( size*sizeof(REAL) );
	u_mean		= (REAL *) malloc ( size*sizeof(REAL) );
	v_mean		= (REAL *) malloc ( size*sizeof(REAL) );
  w_mean		= (REAL *) malloc ( size*sizeof(REAL) );
  temp	    = (REAL *) malloc ( size*sizeof(REAL) );
	temp_s   	= (REAL *) malloc ( size*sizeof(REAL) );
	temp_mean	= (REAL *) malloc ( size*sizeof(REAL) );
	dens		  = (REAL *) malloc ( size*sizeof(REAL) );	
	dens_s  	= (REAL *) malloc ( size*sizeof(REAL) );
	p       	= (REAL *) malloc ( size*sizeof(REAL) ); 
  tmp1      = (REAL *) malloc ( size*sizeof(REAL) );  
  tmp2      = (REAL *) malloc ( size*sizeof(REAL) );  
  tmp3      = (REAL *) malloc ( size*sizeof(REAL) );  
  ap        = (REAL *) malloc ( size*sizeof(REAL) );
  an        = (REAL *) malloc ( size*sizeof(REAL) );
  as        = (REAL *) malloc ( size*sizeof(REAL) );
  aw        = (REAL *) malloc ( size*sizeof(REAL) );
  ae        = (REAL *) malloc ( size*sizeof(REAL) );
  ab        = (REAL *) malloc ( size*sizeof(REAL) );
  af        = (REAL *) malloc ( size*sizeof(REAL) );
  b        = (REAL *) malloc ( size*sizeof(REAL) );
  gx         = (REAL *) malloc ( size*sizeof(REAL) );  
  gy       = (REAL *) malloc ( size*sizeof(REAL) );
  gz       = (REAL *) malloc ( size*sizeof(REAL) );
  ap0      = (REAL *) malloc ( size*sizeof(REAL) );
  pp      = (REAL *) malloc ( size*sizeof(REAL) );
  flagp      = (REAL *) malloc ( size*sizeof(REAL) );
  flagu      = (REAL *) malloc ( size*sizeof(REAL) );
  flagv      = (REAL *) malloc ( size*sizeof(REAL) );
  flagw      = (REAL *) malloc ( size*sizeof(REAL) );
  var	      = (REAL **) malloc ( 38*sizeof(REAL) );
  
  var[X]     = x;
  var[Y]     = y;
  var[Z]     = z;
  var[VX]    = u;
  var[VY]    = v;
  var[VZ]    = w;
  var[VXS]   = u_s;
  var[VYS]   = v_s;
  var[VZS]   = w_s;
  var[VXM]   = u_mean;
  var[VYM]   = v_mean;
  var[VZM]   = w_mean;
  var[DEN]   = dens;
  var[DENS]  = dens_s;
  var[IP]    = p;
  var[TEMP]  = temp;
  var[TEMPS] = temp_s;
  var[TEMPM] = temp_mean;
  var[AP]    = ap;
  var[AN]    = an;
  var[AS]    = as;
  var[AW]    = aw;
  var[AE]    = ae;
  var[AB]    = ab;
  var[AF]    = af;
  var[B]     = b;
  var[TMP1]  = tmp1;
  var[TMP2]  = tmp2;
  var[TMP3]  = tmp3;
  var[GX]    = gx;
  var[GY]    = gy;
  var[GZ]    = gz;
  var[AP0]   = ap0;
  var[PP]    = pp;
  var[FLAGP] =flagp;
  var[FLAGU] =flagu;
  var[FLAGV] =flagv;
  var[FLAGW] =flagw;


	if ( !x || !y || !z || !u || !v || !w || !u_s || !v_s || !w_s || 
       !u_mean || !v_mean || !w_mean || 
       !dens || !dens_s || !temp || !temp_s || !temp_mean || 
       !tmp1 || !tmp2 || !tmp3 ||
       !ap || !ae || !aw || !as || !an || !ab || !af || !b || !gx || !gy || !gz || !ap0 || !pp ) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}

	return ( 1 );
} /** allocate_data() **/


/******************************************************************************
   main --- main routine
******************************************************************************/
int main()
{ 
  int imax, jmax, kmax;
  
  para.geom = &geom;
  para.outp = &outp1;
  para.prob = &prob;
  para.mytime = &mytime;
  para.bc     = &bc;
  para.solv   = &solv;

  initial(&para);

  imax = geom.imax, jmax = geom.jmax, kmax = geom.kmax;

  if(!allocate_data( ))
    exit ( 1 );
   
   clear_data(&para, var);	
   set_grid(&para, var); 
   set_bnd_vel(&para, var, VY, v);

   write_data1(&para, var, "grid");

   FFD_solver(&para, var);


   free_data(var);

   getchar();
  
	exit ( 0 );
} // End of main( )