#include <math.h>

#include "data_structure.h"
#include "solver_gs.h"
#include "boundary.h"

/******************************************************************************
| Gauss-Seidel Solver
******************************************************************************/ 
void GS_P(PARA_DATA *para, REAL **var, int Type, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;

  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;
 
  while(residual>0.000001 && it <5)

   //for(it=0; it<10; it++)
  {
         tmp1=0;
		 tmp2=0.0000000001;
		 it +=1;

    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
		{
		  if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}
    
    for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++)
		{
		  if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}

    for(i=imax; i>=1; i--)
      for(j=jmax; j>=1; j--)
        for(k=1; k<=kmax; k++)
		{
		  if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}

    for(j=jmax; j>=1; j--)
      for(i=imax; i>=1; i--)
        for(k=1; k<=kmax; k++)
		{
		  if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}


   FOR_EACH_CELL
          if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;
    tmp1 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)] 
          - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
          - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
          - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
          - b[IX(i,j,k)]);
    tmp2 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)]);
   END_FOR

	     residual  = tmp1 /tmp2;
   
  }

} // End of Gauss-Seidel( )

void Gauss_Seidel(PARA_DATA *para, REAL **var, int Type, REAL *x)
{
  REAL *as = var[AS], *aw = var[AW], *ae = var[AE], *an = var[AN];
  REAL *ap = var[AP], *af = var[AF], *ab = var[AB], *b = var[B];
  int imax = para->geom->imax, jmax= para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);  
  int i, j, k, it=0;
  float tmp1,tmp2,residual=1;

  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;

  switch (Type)
  {
	  case VX:
		imax=para->geom->imax-1;
		i1= i1-1;
		break;
	  case VY:
        jmax=para->geom->jmax-1;
		j1= j1-1;
		break;
	  case VZ:
		kmax=para->geom->kmax-1;
		break;
  }

  //for(it=0; it<10; it++)
  while(residual>0.000001 && it <5)
  {
         tmp1=0;
		 tmp2=0.0000000001;
		 it += 1;


    for(i=1; i<=imax; i++)
      for(j=1; j<=jmax; j++)
        for(k=1; k<=kmax; k++)
		{

          if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}    
    for(j=1; j<=jmax; j++)
      for(i=1; i<=imax; i++)
        for(k=1; k<=kmax; k++)
		{

          if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}    
    for(i=imax; i>=1; i--)
      for(j=jmax; j>=1; j--)
        for(k=1; k<=kmax; k++)
		{

          if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}    

    for(j=jmax; j>=1; j--)
      for(i=imax; i>=1; i--)
        for(k=1; k<=kmax; k++)
		{

          if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

          x[IX(i,j,k)] = (  ae[IX(i,j,k)]*x[IX(i+1,j,k)] 
                          + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
                          + an[IX(i,j,k)]*x[IX(i,j+1,k)]
                          + as[IX(i,j,k)]*x[IX(i,j-1,k)]
                          + af[IX(i,j,k)]*x[IX(i,j,k+1)]
                          + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
                          + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
		}    


   FOR_EACH_CELL
          if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;
    tmp1 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)] 
          - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
          - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
          - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
          - b[IX(i,j,k)]);
   tmp2 += fabs(ap[IX(i,j,k)]*x[IX(i,j,k)]);
   END_FOR

	     residual  = tmp1 /tmp2;
   
  }

  //printf("*************the residual is %f \t********** \n", residual);

} // End of Gauss-Seidel( )
