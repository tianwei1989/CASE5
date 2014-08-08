/* ----------------------------------------------------------------------------
  
  Filename:	    write.c

  Written by:   Wangda Zuo

	Task:	        output the data as format for tecplot	

---------------------------------------------------------------------------- */	

#include <stdio.h>
#include <string.h>

#include "data_structure.h"
#include "write_data.h"

#include "utility.h"

FILE *file1;

/******************************************************************************
| Write the data to a file for Tecplot 
******************************************************************************/
int write_data(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];
  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 

  int i1 = para->geom->i1, i2 = para->geom->i2; 

  int j1 = para->geom->j1, j2 = para->geom->j2;


  char filename[20];  

  /*change k to degC */
  TempRever(para, var);
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }


    for(j=0; j<=jmax+1; j++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                u[IX(imax+1,j,k)] = u[IX(imax,j,k)];
				um[IX(imax+1,j,k)] = um[IX(imax,j,k)];   

                     for(i=imax; i>=1; i--)
                     {
                          u[IX(i,j,k)] = 0.5f * (u[IX(i,j,k)]+u[IX(i-1,j,k)]);
                          um[IX(i,j,k)] = 0.5f * (um[IX(i,j,k)]+um[IX(i-1,j,k)]);
                      }
            }
	  }


	
    for(i=0; i<=imax+1; i++)
	  {
		for(k=0; k<=kmax+1; k++)
          {
            
                v[IX(i,jmax+1,k)] = v[IX(i,jmax,k)];
                vm[IX(i,jmax+1,k)] = vm[IX(i,jmax,k)];  

                     for(j=jmax; j>=1; j--)
                     {
                          v[IX(i,j,k)] = 0.5f * (v[IX(i,j,k)]+v[IX(i,j-1,k)]);
                          vm[IX(i,j,k)] = 0.5f * (vm[IX(i,j,k)]+vm[IX(i,j-1,k)]);
                      }
            }
	  }
   

	 for(i=0; i<=imax+1; i++)
	  {
		for(j=0; j<=jmax+1; j++)
          {
            
                w[IX(i,j,kmax+1)] = w[IX(i,j,kmax)];
                wm[IX(i,j,kmax+1)] = wm[IX(i,j,kmax)];  

                     for(k=kmax; k>=1; k--)
                     {
                          w[IX(i,j,k)] = 0.5f * (w[IX(i,j,k)]+w[IX(i,j,k-1)]);
                          wm[IX(i,j,k)] = 0.5f * (wm[IX(i,j,k)]+wm[IX(i,j,k-1)]);
                      }
            }
	  }
  
  //FOR_EACH_CELL
	//  if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) T[IX(i,j,k)]=308.3;
  //END_FOR


 //W-S-B
  p[IX(0,0,0)] = (p[IX(0,1,0)]+p[IX(1,0,0)]+p[IX(0,0,1)]) / 3.0f;
  //W-N-B
  p[IX(0,jmax+1,0)] = ( p[IX(1,jmax+1,0)]+p[IX(0,jmax,0)]
                         +p[IX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  p[IX(imax+1,0,0)] = ( p[IX(imax,0,0)]+p[IX(imax+1,1,0)]
                         +p[IX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  p[IX(imax+1,jmax+1,0)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                              +p[IX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  p[IX(0,0,kmax+1)] = ( p[IX(0,1,kmax+1)]+p[IX(1,0,kmax+1)]
                         +p[IX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  p[IX(0,jmax+1,kmax+1)] = ( p[IX(1,jmax+1,kmax+1)]+p[IX(0,jmax,kmax+1)]
                              +p[IX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  p[IX(imax+1,0,kmax+1)] = ( p[IX(imax,0,kmax+1)]+p[IX(imax+1,1,kmax+1)]
                              +p[IX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  p[IX(imax+1,jmax+1,kmax+1)] = ( p[IX(imax,jmax+1,0)]+p[IX(imax+1,jmax,0)]
                                   +p[IX(imax+1,jmax+1,kmax)]) / 3.0f;

								    //W-S-B
  T[IX(0,0,0)] = (T[IX(0,1,0)]+T[IX(1,0,0)]+T[IX(0,0,1)]) / 3.0f;
  //W-N-B
  T[IX(0,jmax+1,0)] = ( T[IX(1,jmax+1,0)]+T[IX(0,jmax,0)]
                         +T[IX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  T[IX(imax+1,0,0)] = ( T[IX(imax,0,0)]+T[IX(imax+1,1,0)]
                         +T[IX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  T[IX(imax+1,jmax+1,0)] = ( T[IX(imax,jmax+1,0)]+T[IX(imax+1,jmax,0)]
                              +T[IX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  T[IX(0,0,kmax+1)] = ( T[IX(0,1,kmax+1)]+T[IX(1,0,kmax+1)]
                         +T[IX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  T[IX(0,jmax+1,kmax+1)] = ( T[IX(1,jmax+1,kmax+1)]+T[IX(0,jmax,kmax+1)]
                              +T[IX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  T[IX(imax+1,0,kmax+1)] = ( T[IX(imax,0,kmax+1)]+T[IX(imax+1,1,kmax+1)]
                              +T[IX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  T[IX(imax+1,jmax+1,kmax+1)] = ( T[IX(imax,jmax+1,0)]+T[IX(imax+1,jmax,0)]
                                   +T[IX(imax+1,jmax+1,kmax)]) / 3.0f;

  
  fprintf( file1, "TITLE = ");

  fprintf( file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, Nx=%d, Ny=%d, Nz=%d \"\n",
           para->mytime->dt, para->mytime->t, para->prob->nu, para->geom->Lx, para->geom->Ly, para->geom->Lz,
           imax+2, jmax+2, kmax+2);

  fprintf( file1, 
           "VARIABLES =X, Y, Z, I, J, K, U, V, W, T, DEN, IP \n");
  fprintf( file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n", imax+2, jmax+2, kmax+2 );

  FOR_ALL_CELL
    fprintf( file1, "%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",
       x[IX(i,j,k)], y[IX(i,j,k)], z[IX(i,j,k)], i, j, k, u[IX(i,j,k)], v[IX(i,j,k)], w[IX(i,j,k)], T[IX(i,j,k)],
       d[IX(i,j,k)], p[IX(i,j,k)]);    
  END_FOR

  fclose(file1);
  
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()

/******************************************************************************
| Caclulate the variabels at 8 corners 
******************************************************************************/
void corners(PARA_DATA *para, REAL **var, REAL *psi)
{
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  //W-S-B
  psi[IX(0,0,0)] = (psi[IX(0,1,0)]+psi[IX(1,0,0)]+psi[IX(0,0,1)]) / 3.0f;
  //W-N-B
  psi[IX(0,jmax+1,0)] = ( psi[IX(1,jmax+1,0)]+psi[IX(0,jmax,0)]
                         +psi[IX(0,jmax+1,1)]) / 3.0f;
  //E-S-B
  psi[IX(imax+1,0,0)] = ( psi[IX(imax,0,0)]+psi[IX(imax+1,1,0)]
                         +psi[IX(imax+1,0,1)]) / 3.0f;
  //E-N-B
  psi[IX(imax+1,jmax+1,0)] = ( psi[IX(imax,jmax+1,0)]+psi[IX(imax+1,jmax,0)]
                              +psi[IX(imax+1,jmax+1,1)]) / 3.0f;
  //W-S-F
  psi[IX(0,0,kmax+1)] = ( psi[IX(0,1,kmax+1)]+psi[IX(1,0,kmax+1)]
                         +psi[IX(0,0,kmax)]) / 3.0f;  
  //W-N-F
  psi[IX(0,jmax+1,kmax+1)] = ( psi[IX(1,jmax+1,kmax+1)]+psi[IX(0,jmax,kmax+1)]
                              +psi[IX(0,jmax+1,kmax)]) / 3.0f;

  //E-S-F
  psi[IX(imax+1,0,kmax+1)] = ( psi[IX(imax,0,kmax+1)]+psi[IX(imax+1,1,kmax+1)]
                              +psi[IX(imax+1,0,kmax)]) / 3.0f;
  //E-N-F
  psi[IX(imax+1,jmax+1,kmax+1)] = ( psi[IX(imax,jmax+1,0)]+psi[IX(imax+1,jmax,0)]
                                   +psi[IX(imax+1,jmax+1,kmax)]) / 3.0f;

} //corners()

int write_data1(PARA_DATA *para, REAL **var, char *name)
{
  int i, j, k;
  int imax=para->geom->imax, jmax=para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int n = para->mytime->t_step;
  int i1=para->geom->i1,i2=para->geom->i2;
  int j1=para->geom->j1,j2=para->geom->j2;
  int k2=para->geom->k2;
  REAL *x = var[X], *y = var[Y], *z =var[Z];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ], *p = var[IP];
  REAL *um = var[VXM], *vm = var[VYM], *wm = var[VZM], *d = var[DEN];
  REAL *T = var[TEMP], *Tm = var[TEMPM];
  REAL *tmp1 = var[TMP1], *tmp2 = var[TMP2], *tmp3 = var[TMP3];

  char filename[20];  
  
  strcpy(filename, name);
  strcat(filename, ".plt");

  /* open output file */
  if((file1 = fopen( filename, "w" ))==NULL)
  {
    fprintf(stderr,"Error:can not open input file!\n");
    return -1;
  }

     
  fprintf( file1, "TITLE = ");

  fprintf( file1, "\"dt=%fs, t=%fs, nu=%f, Lx=%d, Ly=%d, Lz%d, Nx=%d, Ny=%d, Nz=%d \"\n",
           para->mytime->dt, para->mytime->t, para->prob->nu, para->geom->Lx, para->geom->Ly, para->geom->Lz,
           imax+2, jmax+2, kmax+2);

  fprintf( file1, 
           "VARIABLES =X, Y, Z, I, J, K, U, V, W \n");
  fprintf( file1, "ZONE F=POINT, I=%d, J=%d, K=%d\n",  imax+2, jmax+2, kmax+2); //imax+2, jmax+2, kmax+2

  FOR_ALL_CELL
	if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2)
	 continue;
    fprintf( file1, "%f\t%f\t%f\t%d\t%d\t%d\t%f\t%f\t%f\n",
	x[IX(i,j,k)], y[IX(i,j,k)], z[IX(i,j,k)], i, j, k,var[TMP1][IX(i,j,k)],var[TMP2][IX(i,j,k)],var[TMP3][IX(i,j,k)]);    
  END_FOR

  fclose(file1);
  
  printf("The data file %s has been written!\n", name);
  return 0;

} //write_data()