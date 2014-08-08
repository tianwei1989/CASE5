#include "data_structure.h"
#include "utility.h"

/******************************************************************************
| check the residual of equation
******************************************************************************/
REAL check_residual(PARA_DATA *para, REAL **var, REAL *x)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *ap = var[AP], *ab = var[AB], *af = var[AF], *b = var[B];  
  REAL tmp, residual = 0.0; 

  FOR_EACH_CELL
    tmp = ap[IX(i,j,k)]*x[IX(i,j,k)] 
        - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
        - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
        - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
        - b[IX(i,j,k)];
    residual += tmp * tmp;
  END_FOR
    
  return residual / (imax*jmax*kmax);

}// End of check_residual( )

void swap(PARA_DATA *para, REAL **var)
{
	int i, size = (para->geom->imax + 2)*(para->geom->jmax+2)*(para->geom->kmax+2);
  
  
 	for(i=0; i<size; i++) 
  {
    var[VX][i]     = var[TMP1][i];
    var[VY][i]     = var[TMP2][i];
    var[VZ][i]     = var[TMP3][i];
  }
}

/******************************************************************************
| Reverse the k to degC
******************************************************************************/
void TempRever(PARA_DATA *para, REAL **var)
{
  int imax = para->geom->imax, jmax = para->geom->jmax; 
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *tempC = var[TEMP];
  FOR_ALL_CELL
    tempC[IX(i,j,k)] -= 273.15;
  END_FOR
}