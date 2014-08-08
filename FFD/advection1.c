#include <math.h>

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
  REAL x0, y0, z0, x_1, y_1,z_1;
  REAL dt = para->mytime->dt; 
  REAL u0, v0, w0, x_min, x_max, y_min, y_max, z_min, z_max;
  REAL *x = var[X],  *y = var[Y],  *z = var[Z]; 
  REAL *gx = var[GX],  *gy = var[GY],  *gz = var[GZ]; 
  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;

  REAL *u = var[VX], *v = var[VY], *w = var[VZ];

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
    
    x_min = x[IX(0,j,k)];      x_max = x[IX(imax+1,j,k)];  
    //y_min = y[IX(i,0,k)];      y_max = y[IX(i,jmax+1,k)];
    //z_min = z[IX(i,j,0)];      z_max = z[IX(i,j,kmax+1)];
	y_min = y[IX(i,1,k)];      y_max = y[IX(i,jmax,k)];
    z_min = z[IX(i,j,1)];      z_max = z[IX(i,j,kmax)];
        
    x0 =gx[IX(i,j,k)] - u0*dt; 
    y0 = y[IX(i,j,k)] - v0*dt;
    z0 = z[IX(i,j,k)] - w0*dt;
    x0 = x0<x_min ? x_min : x0;
    x0 = x0>x_max ? x_max : x0;
    y0 = y0<y_min ? y_min : y0;
    y0 = y0>y_max ? y_max : y0;
    z0 = z0<z_min ? z_min : z0;
    z0 = z0>z_max ? z_max : z0;

    p = i; q = j; r = k;  
    

/*	if(u0>0 && v0>0 && w0>0)
	{
	  LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x<=0 ||x<=x2) x0=gx[IX(p,q,r)];
			 if(y<=0 ||y<=y2) y0=y[IX(p,q,r)];
			 if(z<=0 ||z<=z1) z0=z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0<gx[IX(p,q,r)]) p -=1;
		     if(y0<y[IX(p,q,r)])  q -=1;
			 if(z0<z[IX(p,q,r)])  r -=1;
		  }
	  }
	}

	else if(u0<0 && v0>0 && w0>0)
	{
	  LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x>=x1 ||x>=Lx) x0=gx[IX(p,q,r)];
			 if(y<=0 ||y<=y2) y0=y[IX(p,q,r)];
			 if(z<=0 ||z<=z1) z0=z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0>gx[IX(p+1,q,r)]) p +=1;
		     if(y0<y[IX(p,q,r)])  q -=1;
			 if(z0<z[IX(p,q,r)])  r -=1;
		  }
	  }
	}
    
	else if(u0>0 && v0<0 && w0>0)
	{
	  LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x<=0 ||x<=x2) x0=gx[IX(p,q,r)];
			 if(y>=y1 ||y>=Ly) y0=y[IX(p,q,r)];
			 if(z<=0 ||z<=z1) z0=z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0<gx[IX(p,q,r)]) p -=1;
		     if(y0>y[IX(p,q+1,r)])  q +=1;
			 if(z0<z[IX(p,q,r)])  r -=1;
		  }
	  }
	}

	else if(u0>0 && v0>0 && w0<0)
	{
	  LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x<=0 ||x<=x2) x0=gx[IX(p,q,r)];
			 if(y<=0 ||y<=y2) y0=y[IX(p,q,r)];
			 if(z>=Lz)        z0=z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0<gx[IX(p,q,r)]) p -=1;
		     if(y0<y[IX(p,q,r)])  q -=1;
			 if(z0>z[IX(p,q,r+1)])  r +=1;
		  }
	  }
	}
	else if(u0>0 && v0<0 && w0<0)
	{
	  LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x<=0 ||x<=x2) x0=gx[IX(p,q,r)];
		     if(y>=y1 ||y>=Ly) y0=y[IX(p,q,r)];
			 if(z>=Lz)        z0=z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0<gx[IX(p,q,r)]) p -=1;
		     if(y0>y[IX(p,q+1,r)])  q +=1;
			 if(z0>z[IX(p,q,r+1)])  r +=1;
		  }
	  }
	}
	else if(u0<0 && v0>0 && w0<0)
	{
	  LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x>=x1 ||x>=Lx) x0=gx[IX(p,q,r)];
			 if(y<=0 ||y<=y2) y0=y[IX(p,q,r)];
			 if(z>=Lz)        z0=z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0>gx[IX(p+1,q,r)]) p +=1;
		     if(y0<y[IX(p,q,r)])  q -=1;
			 if(z0>z[IX(p,q,r+1)])  r +=1;
		  }
	  }
	}

		else if(u0<0 && v0<0 && w0>0)
	{
	  LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x>=x1 ||x>=Lx) x0=gx[IX(p,q,r)];
			 if(y>=y1 ||y>=Ly) y0=y[IX(p,q,r)];
			 if(z<=0 ||z<=z1) z0=z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0>gx[IX(p+1,q,r)]) p +=1;
		     if(y0>y[IX(p,q+1,r)])  q +=1;
			 if(z0<z[IX(p,q,r)])  r -=1;
		  }
	  }
	}
		else
	{ 
		LOCATION=1;
	  while(LOCATION=1)
	  {
		  if(FLAG[IX(p,q,r)]=1)
		  {
		     if(x>=x1 ||x>=Lx) x0=gx[IX(p,q,r)];
			 if(y>=y1 ||y>=Ly) y0= y[IX(p,q,r)];
			 if(z>=Lz)         z0= z[IX(p,q,r)];
			 LOCATION=0;
	 	  }
		  else
		  {  if(x0>gx[IX(p+1,q,r)]) p +=1;
		     if(y0>y[IX(p,q+1,r)])  q +=1;
			 if(z0>z[IX(p,q,r+1)])  r +=1;
		  }
	  }
			
	}

	*/











    /*-----------------------------------------------------------------------
    | Step 2: Decide Index and Interpolate
    -----------------------------------------------------------------------*/
    if(u0 > 0)
      while(x0<gx[IX(p,j,k)] ) p -= 1;        
    else 
      while(x0>gx[IX(p+1,j,k)] ) p += 1;

    if(v0 > 0)
      while(y0<y[IX(i,q,k)] ) q -= 1;    
    else
      while(y0>y[IX(i,q+1,k)] ) q += 1;

    if(w0 > 0)
      while(z0<z[IX(i,j,r)] ) r -= 1;
    else
      while(z0>z[IX(i,j,r+1)] ) r += 1;
    
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
    
   // x_min = x[IX(0,j,k)];      x_max = x[IX(imax+1 ,j,k)];  
    y_min = y[IX(i,0,k)];      y_max = y[IX(i,jmax+1,k)];
    //z_min = z[IX(i,j,0)];      z_max = z[IX(i,j,kmax+1 )];

	 x_min = x[IX(1,j,k)];      x_max = x[IX(imax ,j,k)];
	 z_min = z[IX(i,j,1)];     z_max = z[IX(i,j,kmax )];
        
    x0 = x[IX(i,j,k)] - u0*dt; 
    y0 = gy[IX(i,j,k)] - v0*dt;
    z0 = z[IX(i,j,k)] - w0*dt;
    x0 = x0<x_min ? x_min : x0;
    x0 = x0>x_max ? x_max : x0;
    y0 = y0<y_min ? y_min : y0;
    y0 = y0>y_max ? y_max : y0;
    z0 = z0<z_min ? z_min : z0;
    z0 = z0>z_max ? z_max : z0;

    p = i; q = j; r = k;  
    /*-----------------------------------------------------------------------
    | Step 2: Decide Index and Interpolate
    -----------------------------------------------------------------------*/
    if(u0 > 0)
      while(x0<x[IX(p,j,k)] ) p -= 1;        
    else 
      while(x0>x[IX(p+1,j,k)] ) p += 1;

    if(v0 > 0)
      while(y0<gy[IX(i,q,k)] ) q -= 1;    
    else
      while(y0>gy[IX(i,q+1,k)]) q += 1;

    if(w0 > 0)
      while(z0<z[IX(i,j,r)] ) r -= 1;
    else
      while(z0>z[IX(i,j,r+1)] ) r += 1;
    
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
    
   // x_min = x[IX(0,j,k)];      x_max = x[IX(imax+1 ,j,k)];  
   // y_min = y[IX(i,0,k)];      y_max = y[IX(i,jmax+1 ,k)];
    z_min = z[IX(i,j,0)];      z_max = z[IX(i,j,kmax+1)];

	x_min = x[IX(1,j,k)];      x_max = x[IX(imax ,j,k)];  
    y_min = y[IX(i,1,k)];      y_max = y[IX(i,jmax ,k)];
        
    x0 = x[IX(i,j,k)] - u0*dt; 
    y0 = y[IX(i,j,k)] - v0*dt;
    z0 = gz[IX(i,j,k)] - w0*dt;
    x0 = x0<x_min ? x_min : x0;
    x0 = x0>x_max ? x_max : x0;
    y0 = y0<y_min ? y_min : y0;
    y0 = y0>y_max ? y_max : y0;
    z0 = z0<z_min ? z_min : z0;
    z0 = z0>z_max ? z_max : z0;

    p = i; q = j; r = k;  
    /*-----------------------------------------------------------------------
    | Step 2: Decide Index and Interpolate
    -----------------------------------------------------------------------*/
    if(u0 > 0)
      while(x0<x[IX(p,j,k)] ) p -= 1;        
    else 
      while(x0>x[IX(p+1,j,k)] ) p += 1;

    if(v0 > 0)
      while(y0<y[IX(i,q,k)] ) q -= 1;    
    else
      while(y0>y[IX(i,q+1,k)] ) q += 1;

    if(w0 > 0)
      while(z0<gz[IX(i,j,r)]) r -= 1;
    else
      while(z0>gz[IX(i,j,r+1)] ) r += 1;
    
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


