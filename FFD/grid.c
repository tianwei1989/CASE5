#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "data_structure.h"
#include "grid.h"

/******************************************************************************
| Set Grid
******************************************************************************/
void set_grid(PARA_DATA *para, REAL **var)
{

  if(para->geom->uniform)
    set_uniform_grid(para, var);
  else 
    set_nonuniform_grid(para,var);
} // End of set_grid( )

/******************************************************************************
| Set Uniform Grid
******************************************************************************/
void set_uniform_grid(PARA_DATA *para, REAL **var)
{
  int i, j, k;  
  int imax = para->geom->imax, jmax = para->geom->jmax;  
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  REAL dx1,dx2,dx3; 
  REAL dy1,dy2,dy3;
  REAL dz1,dz2,dz3,dz4;
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 
  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  REAL z1 = para->geom->z1, z2 = para->geom->z2 ,z3 = para->geom->z3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  REAL x1 = para->geom->x1, x2 = para->geom->x2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;
  REAL y1 = para->geom->y1, y2 = para->geom->y2;
  
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];

  switch(para->solv->caseID)
  { 
  case 1:


   dx1=Lx/(REAL) imax;
   dy1=Ly/(REAL) jmax;
   dz1=z1/ (REAL) k1;
   dz2= (z2-z1)/(REAL) (k2-k1);
   dz3= (Lz-z2)/(REAL) (kmax-k2);
 
    /*---------------------------------------------------------------------------
  | Central Domain
  ---------------------------------------------------------------------------*/
  FOR_ALL_CELL 

   gx[IX(i,j,k)] = dx1 * (REAL) i;
   gy[IX(i,j,k)] = dy1 * (REAL) j;
  
   if(k<=k1)  gz[IX(i,j,k)] = dz1 * (REAL) k;
   else if (k<=k2)   gz[IX(i,j,k)] = z1 + dz2 * (REAL) (k-k1);
   else  gz[IX(i,j,k)] = z2 + dz3 * (REAL) (k-k2);
     
	  if(i<1)  x[IX(i,j,k)]= 0;
	  else if (i>imax) x[IX(i,j,k)]= Lx;
	  else    x[IX(i,j,k)]= 0.5f* (gx[IX(i,j,k)]+gx[IX(i-1,j,k)]);
	 
	  if(j<1)  y[IX(i,j,k)]= 0;
	  else if(j>jmax) y[IX(i,j,k)]= Ly;
	  else y[IX(i,j,k)]= 0.5f* (gy[IX(i,j,k)]+gy[IX(i,j-1,k)]);

	  if(k<1)  z[IX(i,j,k)]= 0;
	  else if(k>kmax) z[IX(i,j,k)]= Lz;
	  else z[IX(i,j,k)]= 0.5f* (gz[IX(i,j,k)]+gz[IX(i,j,k-1)]);
 	
   END_FOR    

		   break;


  case 2:

   dx1=x1/(REAL) i1;
   dx2=(x2-x1)/(REAL) (i2-i1);
   dx3=(Lx-x2)/(REAL) (imax-i2);
   dy1=y1/(REAL) j1;
   dy2=(y2-y1)/(REAL) (j2-j1);
   dy3=(Ly-y2)/(REAL) (jmax-j2);
   dz1= z1/(REAL) k1;
   dz2=(z2-z1)/(REAL) (k2-k1);
   dz3=(z3-z2)/(REAL) (k3-k2);
   dz4=(Lz-z3)/(REAL) (kmax-k3);
   
    /*---------------------------------------------------------------------------
  | Central Domain
  ---------------------------------------------------------------------------*/
  FOR_ALL_CELL   

	
   if(i<=i1)  gx[IX(i,j,k)] = dx1 * (REAL) i;
   else if (i<=i2)   gx[IX(i,j,k)] = x1 + dx2 * (REAL) (i-i1);
   else  gx[IX(i,j,k)] = x2 + dx3 * (REAL) (i-i2);

   if(j<=j1)  gy[IX(i,j,k)] = dy1 * (REAL) j;
   else if (j<=j2)   gy[IX(i,j,k)] = y1 + dy2 * (REAL) (j-j1);
   else  gy[IX(i,j,k)] = y2 + dy3 * (REAL) (j-j2);


   if(k<=k1)  gz[IX(i,j,k)] = dz1 * (REAL) k;
   else if (k<=k2)   gz[IX(i,j,k)] = z1 + dz2 * (REAL) (k-k1);
   else if (k<=k3)   gz[IX(i,j,k)] = z2 + dz3 * (REAL) (k-k2); 
   else     gz[IX(i,j,k)] = z3 + dz4 * (REAL) (k-k3);

     
	  if(i<1)  x[IX(i,j,k)]= 0;
	  else if (i>imax) x[IX(i,j,k)]= Lx;
	  else    x[IX(i,j,k)]= 0.5f* (gx[IX(i,j,k)]+gx[IX(i-1,j,k)]);
	 
	  if(j<1)  y[IX(i,j,k)]= 0;
	  else if(j>jmax) y[IX(i,j,k)]= Ly;
	  else y[IX(i,j,k)]= 0.5f* (gy[IX(i,j,k)]+gy[IX(i,j-1,k)]);

	  if(k<1)  z[IX(i,j,k)]= 0;
	  else if(k>kmax) z[IX(i,j,k)]= Lz;
	  else z[IX(i,j,k)]= 0.5f* (gz[IX(i,j,k)]+gz[IX(i,j,k-1)]);


   if(i<1 || j<1 || k<1 || i>imax || j>jmax || k>kmax) flagp[IX(i,j,k)]=1;
   if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) flagp[IX(i,j,k)]=1;
   
   if(i<1 || j<1 || k<1 || i>imax-1 || j>jmax ||k>kmax) flagu[IX(i,j,k)]=1;
   if(i>i1-1 && i<=i2 && j>j1 && j<=j2 && k<=k2) flagu[IX(i,j,k)]=1;

   if(i<1 || j<1 || k<1 || i>imax || j>jmax-1 ||k>kmax) flagv[IX(i,j,k)]=1;
   if(i>i1 && i<=i2 && j>j1-1 && j<=j2 && k<=k2) flagv[IX(i,j,k)]=1;

   if(i<1 || j<1 || k<1 || i>imax || j>jmax ||k>kmax-1) flagw[IX(i,j,k)]=1;
   if(i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) flagw[IX(i,j,k)]=1;

   if(j<1 && k>k3) {flagv[IX(i,j,k)]=0; flagp[IX(i,j,k)]=0;}

 
  END_FOR  

    
	  break;
  }
          
} // End of set_uniform_grid( )

/******************************************************************************
| Set Nonuniform Grid
******************************************************************************/
void set_nonuniform_grid(PARA_DATA *para, REAL **var)
{
  int i, j, k;  
  int imax = para->geom->imax, jmax = para->geom->jmax;    
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2); 
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz;
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  int  k2= para->geom->k2, k1 = para->geom->k1;
  REAL z2 = para->geom->z2, z1 = para->geom->z1; 
  REAL fac, x1, y_temp;
   
      fac = 1.2f;
     
      for(i=0; i<=imax; i++)
      {        
        x1 = Lx * (0.5 + 0.5 * tanh(fac*(2.0*i/imax-1.0))/tanh(fac));     
        for(j=0; j<=jmax+1; j++) 
          for(k=0; k<=kmax+1; k++) 
            gx[IX(i,j,k)] = x1;     
      }  


	  fac = 1.2f;
     
      for(j=0; j<=jmax; j++)
      {        
        x1 = Ly * (0.5 + 0.5 * tanh(fac*(2.0*j/jmax-1.0))/tanh(fac));     
        for(i=0; i<=imax+1; i++) 
          for(k=0; k<=kmax+1; k++) 
            gy[IX(i,j,k)] = x1;     
      }  
      
    	 

	fac = 1.1f;

	for(i=0; i<=imax+1; i++)  
	{
        for(j=0; j<=jmax+1; j++)
        {
               
          for(k=0; k<=k1; k++)  gz[IX(i,j,k)] = z1 / (REAL) k1 * (REAL)k;
                  
          for(k=k1+1; k<=k2; k++) gz[IX(i,j,k)] = z1 + (z2-z1)* (0.5f + 0.5f * tanh(fac*(2.0f* (REAL)(k-k1)/ (REAL)(k2-k1)-1.0))/tanh(fac));

		  for(k=k2+1;k<=kmax;k++) gz[IX(i,j,k)] =z2+ ( Lz-z2)/ (REAL) (kmax-k2) * (REAL)(k-k2);
        }
	}

	  FOR_ALL_CELL

	  if(i<1)  x[IX(i,j,k)]= 0;
	  else if (i>imax)  x[IX(i,j,k)]= Lx;
	  else  x[IX(i,j,k)]= 0.5f* (gx[IX(i,j,k)]+gx[IX(i-1,j,k)]);  
	 
	  if(j<1)  y[IX(i,j,k)]= 0;
	  else if(j>jmax) y[IX(i,j,k)]= Ly;
	  else y[IX(i,j,k)]= 0.5f* (gy[IX(i,j,k)]+gy[IX(i,j-1,k)]);

	  if(k<1)  z[IX(i,j,k)]= 0;
	  else if(k>kmax) z[IX(i,j,k)]= Lz;
	  else z[IX(i,j,k)]= 0.5f* (gz[IX(i,j,k)]+gz[IX(i,j,k-1)]);

	  END_FOR

      

} // End of set_nonuniform_grid()