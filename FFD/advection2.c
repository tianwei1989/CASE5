#include <math.h>
#include <stdio.h>
#include "data_structure.h"
#include "boundary.h"
#include "advection.h"
#include "solver.h"
#include "utility.h"
#include "interpolation.h"

/******************************************************************************
| advection
******************************************************************************/
void advection(PARA_DATA *para, REAL **var, int var_type, REAL *d, REAL *d0)
{
	int i, j, k, p, q, r;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int LOCATION;
  REAL x0, y0, z0, x_1, y_1,z_1;
  REAL dt = para->mytime->dt; 
  REAL u0, v0, w0, x_min, x_max, y_min, y_max, z_min, z_max;
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 

  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *flagp = var[FLAGP],*flagu = var[FLAGU],*flagv = var[FLAGV],*flagw = var[FLAGW];
  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;
  REAL z1 = para->geom->z1, z2 = para->geom->z2 ,z3 = para->geom->z3; 
  REAL x1 = para->geom->x1, x2 = para->geom->x2; 
  REAL y1 = para->geom->y1, y2 = para->geom->y2;
  REAL Lx = para->geom->Lx, Ly = para->geom->Ly, Lz = para->geom->Lz; 


  switch (var_type)
  {
	  
  case VX:
   {
  /*---------------------------------------------------------------------------
  | Tracing Back and Interplotaing
  ---------------------------------------------------------------------------*/
  FOR_U_CELL

	if (i>=i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;
    /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    u0 = u[IX(i,j,k)];

    v0 = 0.5f *((v[IX(i,  j,k)]+ v[IX(i,  j-1,k)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
		       +(v[IX(i+1,j,k)]+ v[IX(i+1,j-1,k)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)])) /(x[IX(i+1,j,k)]-x[IX(i,j,k)]);
	
	w0 = 0.5f *((w[IX(i,  j,k)]+ w[IX(i  ,j, k-1)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
		       +(w[IX(i+1,j,k)]+ w[IX(i+1,j, k-1)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)]))/(x[IX(i+1,j,k)]-x[IX(i,j,k)]); 
    
        
    x0 =gx[IX(i,j,k)] - u0*dt; 
    y0 = y[IX(i,j,k)] - v0*dt;
    z0 = z[IX(i,j,k)] - w0*dt;
   
	
	p = i; q = j; r = k;  

  if(u0>0)
  {
	  while(x0<gx[IX(p,q,r)] && flagu[IX(p,q,r)]<0)  p -=1;
      if(flagu[IX(p,q,r)]>0 && x0< gx[IX(p,q,r)]) x0=gx[IX(p,q,r)];
  }
  else
  {
     while(x0>gx[IX(p+1,q,r)] && flagu[IX(p+1,q,r)]<0)  p +=1;
	 if(flagu[IX(p+1,q,r)]>0 && x0>gx[IX(p+1,q,r)])  x0=gx[IX(p+1,q,r)];
  }

  if(v0>0)
  {
	  while(y0<y[IX(p,q,r)] && flagu[IX(p,q,r)]<0)  q -=1;
      if(y0<y[IX(p,q,r)] && flagu[IX(p,q,r)]>0) y0=y[IX(p,q+1,r)];
  }
  else
  {
   while(y0>y[IX(p,q+1,r)] && flagu[IX(p,q+1,r)]<0)  q +=1;
      if(y0>y[IX(p,q+1,r)] && flagu[IX(p,q+1,r)]>0)  y0=y[IX(p,q,r)];
  }

    if(w0>0)
  {
	  while(z0<z[IX(p,q,r)] && flagu[IX(p,q,r)]<0)  r -=1;
      if(z0<z[IX(p,q,r)] && flagu[IX(p,q,r)]>0) z0=z[IX(p,q,r+1)];
  }
  else
  {
   while(z0>z[IX(p,q,r+1)] && flagu[IX(p,q,r+1)]<0)  r +=1;
    if(z0>z[IX(p,q,r+1)] && flagu[IX(p,q,r+1)]>0)  z0=z[IX(p,q,r)];
  }

	     
    x_1 = (x0-gx[IX(p,q,r)]) /(gx[IX(p+1,  q ,r  )]-gx[IX(p,q,r)]); 
    y_1 = (y0- y[IX(p,q,r)]) /(y [IX(p,  q+1, r  )]- y[IX(p,q,r)]);
    z_1 = (z0- z[IX(p,q,r)]) /(z [IX(p,    q, r+1)]- z[IX(p,q,r)]);
             
    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, p, q, r);
	END_FOR

  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d);

   }
	 break;


  case VY:
  {
	/*---------------------------------------------------------------------------
  | Tracing Back and Interplotaing
  ---------------------------------------------------------------------------*/
  FOR_V_CELL

    if (i>i1 && i<=i2 && j>=j1 && j<=j2 && k<=k2) continue;
    /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    
	u0 = 0.5f* ((u[IX(i,j,k)]   + u[IX(i-1,j,  k)]) *(y [IX(i,j+1,k)]-gy[IX(i,j,k)])
		       +(u[IX(i,j+1,k)] + u[IX(i-1,j+1,k)]) *(gy[IX(i,j,  k)]-y[IX(i,j,k)]))/(y[IX(i,j+1,k)]-y[IX(i,j,k)]);
	
    v0 = v[IX(i,j,k)]; 


   	w0 = 0.5f *((w[IX(i,j,  k)]+ w[IX(i,j  ,k-1)])* (y [IX(i,j+1,k)]- gy[IX(i,j,k)])
		       +(w[IX(i,j+1,k)]+ w[IX(i,j+1,k-1)])* (gy[IX(i,j,  k)]- y [IX(i,j,k)]))/(y[IX(i,j+1,k)]-y[IX(i,j,k)]); 
    
       
    x0 = x[IX(i,j,k)] - u0*dt; 
    y0 = gy[IX(i,j,k)] - v0*dt;
    z0 = z[IX(i,j,k)] - w0*dt;
   
    p = i; q = j; r = k;  


	if(u0>0)
  {
	  while(x0<x[IX(p,q,r)] && flagv[IX(p,q,r)]<0)  p -=1;
      if(flagv[IX(p,q,r)]>0 && x0< x[IX(p,q,r)]) x0=x[IX(p+1,q,r)];
  }
  else
  {
     while(x0>x[IX(p+1,q,r)] && flagv[IX(p+1,q,r)]<0)  p +=1;
	 if(flagv[IX(p+1,q,r)]>0 && x0>x[IX(p+1,q,r)])  x0=x[IX(p,q,r)];
  }

  if(v0>0)
  {
	  while(y0<gy[IX(p,q,r)] && flagv[IX(p,q,r)]<0)  q -=1;
      if(y0<gy[IX(p,q,r)] && flagv[IX(p,q,r)]>0) y0=gy[IX(p,q,r)];
  }
  else
  {
   while(y0>gy[IX(p,q+1,r)] && flagv[IX(p,q+1,r)]<0)  q +=1;
      if(y0>gy[IX(p,q+1,r)] && flagv[IX(p,q+1,r)]>0)  y0=gy[IX(p,q+1,r)];
  }

    if(w0>0)
  {
	  while(z0<z[IX(p,q,r)] && flagv[IX(p,q,r)]<0)  r -=1;
      if(z0<z[IX(p,q,r)] && flagv[IX(p,q,r)]>0) z0=z[IX(p,q,r+1)];
  }
  else
  {
   while(z0>z[IX(p,q,r+1)] && flagv[IX(p,q,r+1)]<0)  r +=1;
    if(z0>z[IX(p,q,r+1)] && flagv[IX(p,q,r+1)]>0)  z0=z[IX(p,q,r)];
  }
  

    
    x_1 = (x0- x[IX(p,q,r)])/(x [IX(p+1,q,    r)]-x [IX(p,q,r)]); 
    y_1 = (y0-gy[IX(p,q,r)])/(gy[IX(p,  q+1,  r)]-gy[IX(p,q,r)]);
    z_1 = (z0- z[IX(p,q,r)])/(z [IX(p,  q,  r+1)]-z [IX(p,q,r)]);
             
    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, p, q, r);
	END_FOR

  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d);

		  }
      break;

  case VZ:
   {
		  	/*---------------------------------------------------------------------------
  | Tracing Back and Interplotaing
  ---------------------------------------------------------------------------*/
  FOR_W_CELL

	 if (i>i1 && i<=i2 && j>j1 && j<=j2 && k<=k2) continue;

    /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    
	u0 = 0.5f*((u[IX(i,j,k  )]  + u[IX(i-1,j,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
		      +(u[IX(i,j,k+1)]  + u[IX(i-1,j,k+1)])*(gz[IX(i,j,k  )]- z[IX(i,j,k)]))/(z[IX(i,j,k+1)]-z[IX(i,j,k)]);

	v0 = 0.5f*((v[IX(i,j,k  )]  + v[IX(i,j-1,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
		      +(v[IX(i,j,k+1)]  + v[IX(i,j-1,k+1)])*(gz[IX(i,j,k  )]-z [IX(i,j,k)]))/(z[IX(i,j,k+1)]-z[IX(i,j,k)]); 

	w0 = w[IX(i,j,k)]; 
          
    x0 = x[IX(i,j,k)] - u0*dt; 
    y0 = y[IX(i,j,k)] - v0*dt;
    z0 = gz[IX(i,j,k)] - w0*dt;
   
    p = i; q = j; r = k;  


	if(u0>0)
  {
	  while(x0<x[IX(p,q,r)] && flagw[IX(p,q,r)]<0)  p -=1;
      if(flagw[IX(p,q,r)]>0 && x0< x[IX(p,q,r)]) x0=x[IX(p+1,q,r)];
  }
  else
  {
     while(x0>x[IX(p+1,q,r)] && flagw[IX(p+1,q,r)]<0)  p +=1;
	 if(flagw[IX(p+1,q,r)]>0 && x0>x[IX(p+1,q,r)])  x0=x[IX(p,q,r)];
  }

  if(v0>0)
  {
	  while(y0<y[IX(p,q,r)] && flagw[IX(p,q,r)]<0)  q -=1;
      if(y0<y[IX(p,q,r)] && flagw[IX(p,q,r)]>0) y0=y[IX(p,q+1,r)];
  }
  else
  {
   while(y0>y[IX(p,q+1,r)] && flagw[IX(p,q+1,r)]<0)  q +=1;
      if(y0>y[IX(p,q+1,r)] && flagw[IX(p,q+1,r)]>0)  y0=y[IX(p,q,r)];
  }

    if(w0>0)
  {
	  while(z0<gz[IX(p,q,r)] && flagw[IX(p,q,r)]<0)  r -=1;
      if(z0<gz[IX(p,q,r)] && flagw[IX(p,q,r)]>0) z0=gz[IX(p,q,r)];
  }
  else
  {
   while(z0>gz[IX(p,q,r+1)] && flagw[IX(p,q,r+1)]<0)  r +=1;
    if(z0>gz[IX(p,q,r+1)] && flagw[IX(p,q,r+1)]>0)  z0=gz[IX(p,q,r+1)];
  }

     
    x_1 = (x0- x[IX(p,q,r)])/( x[IX(p+1,q,   r  )]- x[IX(p,q,r)]); 
    y_1 = (y0- y[IX(p,q,r)])/( y[IX(p,  q+1, r  )]- y[IX(p,q,r)]);
    z_1 = (z0-gz[IX(p,q,r)])/(gz[IX(p,  q,   r+1)]-gz[IX(p,q,r)]);
             
    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    d[IX(i,j,k)] = interpolation(para, d0, x_1, y_1, z_1, p, q, r);
	END_FOR
		  
  /*---------------------------------------------------------------------------
  | define the b.c.
  ---------------------------------------------------------------------------*/
  set_bnd(para, var, var_type, d);
		  
	 }
		  break;
  }
		 

} // End of semi_Lagrangian( )


