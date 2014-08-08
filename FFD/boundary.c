#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "data_structure.h"
#include "boundary.h"

#include "inlet_profile.h"

/******************************************************************************
|  Set the boundary conditions
******************************************************************************/
void set_bnd(PARA_DATA *para, REAL **var, int var_type, REAL *psi)
{
  switch(var_type)
  {
    case VX:
      set_bnd_vel(para, var, VX, psi); break;
    case VY:
      set_bnd_vel(para, var, VY, psi); break;
    case VZ:
      set_bnd_vel(para, var, VZ, psi); break;
	case TEMP:
	  set_bnd_temp(para, var, TEMP, psi); break;
 
  }
} // set_bnd() 


/**************************************************************************************************************
|  Set the boundary conditions for velocity
*******************************************************************************************************************/
void set_bnd_vel(PARA_DATA *para, REAL **var, int var_type, REAL *psi)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *x = var[X], *z = var[Z];
  REAL *af = var[AF], *ab = var[AB];
  int caseID = para->solv->caseID;

  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;

  /* --------------------------------------------------------------------------
  | North Boundary
  -------------------------------------------------------------------------- */
  j = jmax; //outflow//
  switch(var_type)
      {
        case VX:
		case VZ:
			  FOR_KI
				  if(k<=k1)
				  {
                  psi[IX(i,j+1,k)] = psi[IX(i,j ,k)];
				  an[IX(i,j,k)] = 0.0;
				  }
				  else  psi[IX(i,j+1,k)] = 0; //nonslip wall
				
              END_FOR
			
			break;
        case VY:
			  FOR_KI
				  if(k<=k1)
				  {
                 // psi[IX(i,j,k)] = psi[IX(i,j-1,k)];
                 // an[IX(i,j-1,k)] = 0.0;
				  psi[IX(i,j,k)] = psi[IX(i,j,k)];
				  //if (psi[IX(i,j,k)]<0)  {psi[IX(i,j,k)]=0;psi[IX(i,j-1,k)]=0;}
				  }
				  else  psi[IX(i,j,k)] =0; //nonslip wall

               END_FOR
			break;
      } 
  // End of North Boundary

 /* --------------------------------------------------------------------------
  | South Boundary
  -------------------------------------------------------------------------- */
  j = 1;   //Inflow//
  switch(var_type)
      {
        case VX:
			   	  FOR_KI
				  if(k>k3)  psi[IX(i,j-1,k)] = 0;
				  else      psi[IX(i,j-1,k)] =0; //nonslip wall
				  END_FOR
			break;
        case VY:
				  FOR_KI
					  if(k>k3)  psi[IX(i,j-1,k)] = para->bc->VY_bcS;
					  else      psi[IX(i,j-1,k)] = 0; //nonslip wall
                  END_FOR
	    	break;
		case VZ: 	
				  FOR_KI
				  if(k>k3)   psi[IX(i,j-1,k)] = 0;
				  else       psi[IX(i,j-1,k)] = 0; //nonslip wall
                  END_FOR
		    break;
      
      } 
          
  // End of South Boundary

  /* --------------------------------------------------------------------------
  | East Boundary
  -------------------------------------------------------------------------- */
  i = imax;  //noslip//
  switch(var_type)
      {
        case VX:
		     	  FOR_JK
                  psi[IX(i,j,k)] = 0;
                  END_FOR
			break;
        case VY:
				  FOR_JK
                   psi[IX(i+1,j,k)] = 0;
				   // psi[IX(i+1,j,k)] = psi[IX(i,j,k)];
					//ae[IX(i,j,k)]=0;
                  END_FOR
			break;
		case VZ: 	
				  FOR_JK
                   psi[IX(i+1,j,k)] = 0;
				  // psi[IX(i+1,j,k)] = psi[IX(i,j,k)];
				  // ae[IX(i,j,k)]=0;
                  END_FOR
			break;
        } 
 
 // End of East Boundary

  /* --------------------------------------------------------------------------
  | West Boundary
  -------------------------------------------------------------------------- */
  i = 1; //nonslip//
  	switch(var_type)
      {
		case VX:
			FOR_JK
            psi[IX(i-1,j,k)] =0;
            END_FOR
		break;
		case VY:
            FOR_JK
            psi[IX(i-1,j,k)] =0;
			//psi[IX(i-1,j,k)] = psi[IX(i,j,k)];
			//aw[IX(i,j,k)]=0;
			END_FOR
		break;
		case VZ:
            FOR_JK
            psi[IX(i-1,j,k)] =0;
			//psi[IX(i-1,j,k)] = psi[IX(i,j,k)];
			//aw[IX(i,j,k)]=0;
			END_FOR
		break;
	   }

 // End of West Boundary

  /* --------------------------------------------------------------------------
  | Back Boundary
  -------------------------------------------------------------------------- */
  k = 1; //noslip//
   
      switch(var_type)
        {
	     case VX:
			      FOR_IJ
                  psi[IX(i,j,k-1)] = 0;
				  END_FOR
			  break;
		 case VY:
			     FOR_IJ
                  psi[IX(i ,j,k-1)] = 0;
				  END_FOR
			 break;
		 case VZ:
			      FOR_IJ
                  psi[IX(i ,j,k-1)] = 0;
				  END_FOR
		     break;
	    }

// End of Back Boundary

  /* --------------------------------------------------------------------------
  | Front Boundary
  -------------------------------------------------------------------------- */
  k = kmax; //noslip//

      switch(var_type)
          {
	     case VX:
			      FOR_IJ
                  psi[IX(i,j,k+1)] = 0;
				  END_FOR
			  break;
		 case VY:
			      FOR_IJ
                  psi[IX(i,j,k+1)] = 0;
				  END_FOR
			 break;
		 case VZ:
			      FOR_IJ
                  psi[IX(i,j,k)] = 0;
				  END_FOR
			 break;
         } 
	
 // End of Front Boundary

  /* --------------------------------------------------------------------------
  | BOX SOUTH
  -------------------------------------------------------------------------- */
  j = j1; //no slip//
  switch(var_type)
      {
        case VX:
			for(i=i1; i<=i2; i++)
				for(k=1;k<=k2; k++)
					//psi[IX(i,j+1,k)] = -psi[IX(i,j,k)];
					psi[IX(i,j+1,k)] = 0;
			break;
		case VZ:
			 for(i=i1+1; i<=i2; i++)
				for(k=1;k<=k2; k++)
				   // psi[IX(i,j+1,k)] = -psi[IX(i,j,k)]; //nonslip wall	
				   psi[IX(i,j+1,k)]  = 0;
			break;
        case VY:
			 for(i=i1+1; i<=i2; i++)
				for(k=1;k<=k2; k++)
				   psi[IX(i,j,k)] =0; //nonslip wall

			break;
      } 
  // End of North Boundary

 /* --------------------------------------------------------------------------
  | BOX NORTH
  -------------------------------------------------------------------------- */
  j = j2;   //no slip
  switch(var_type)
      {
        case VX:
			 for(i=i1; i<=i2; i++)
				for(k=1;k<=k2; k++)
				  //  psi[IX(i,j,k)] = - psi[IX(i,j+1,k)]; //nonslip wall
			         psi[IX(i,j,k)] = 0;
			 break;
        case VY:
			 for(i=i1+1; i<=i2; i++)
				for(k=1;k<=k2; k++)
				    psi[IX(i,j,k)] = 0; //nonslip wall
              
	    	break;
		case VZ: 	
			 for(i=i1+1; i<=i2; i++)
				for(k=1;k<=k2; k++)
			 
			     // psi[IX(i,j,k)] = - psi[IX(i,j+1,k)]; //nonslip wall
			       psi[IX(i,j,k)] = 0;
          
		    break;
      
      } 
          
  // End of South Boundary

  /* --------------------------------------------------------------------------
  | BOX EAST
  -------------------------------------------------------------------------- */
  i = i2;  //noslip//
  switch(var_type)
      {
        case VX:
			 for(j=j1+1;j<=j2; j++)
				for(k=1;k<=k2; k++)
                  psi[IX(i,j,k)] = 0;
              
			break;
        case VY:
			 for(j=j1;j<=j2; j++)
				for(k=1;k<=k2; k++)
                  // psi[IX(i,j,k)] = - psi[IX(i+1,j,k)];
			       psi[IX(i,j,k)] = 0;
			
			break;
		case VZ: 	
			 for(j=j1+1;j<=j2; j++)
				for(k=1;k<=k2; k++)
                 //  psi[IX(i,j,k)] =  - psi[IX(i+1,j,k)];
			        psi[IX(i,j,k)] =  0;
				 
			break;
        } 
 
 // End of East Boundary

  /* --------------------------------------------------------------------------
  | West Boundary
  -------------------------------------------------------------------------- */
  i = i1; //nonslip//
  	switch(var_type)
      {
		case VX:
			 for(j=j1+1;j<=j2; j++)
				for(k=1;k<=k2; k++)
                    psi[IX(i,j,k)] =0;
          
		break;
		case VY:
			 for(j=j1;j<=j2; j++)
				for(k=1;k<=k2; k++)
                  // psi[IX(i,j,k)] = - psi[IX(i-1,j,k)];
			 psi[IX(i,j,k)] = 0;
			
		break;
		case VZ:
             for(j=j1+1;j<=j2; j++)
				for(k=1;k<=k2; k++)
               //     psi[IX(i,j,k)] =- psi[IX(i-1,j,k)];
			  psi[IX(i,j,k)] =0;
	
		break;
	   }

 // End of West Boundary

  /* --------------------------------------------------------------------------
  | BOX TOP
  -------------------------------------------------------------------------- */
  k = k2; //noslip//
   
      switch(var_type)
        {
	     case VX:
			 for(i=i1; i<=i2; i++)
				for(j=j1+1;j<=j2; j++)
                 //  psi[IX(i,j,k)] = -psi[IX(i,j,k+1)];
				 psi[IX(i,j,k)] = 0;
			  break;
		 case VY:
			  for(i=i1+1; i<=i2; i++)
				for(j=j1;j<=j2; j++)
               //   psi[IX(i ,j,k)] =  -psi[IX(i,j,k+1)];
			     psi[IX(i ,j,k)] = 0;
			 break;
		 case VZ:
			  for(i=i1+1; i<=i2; i++)
				for(j=j1+1;j<=j2; j++)
                  psi[IX(i ,j,k)] = 0;
				
		     break;
	    }
 
} // End of set_bnd_vel( )


void set_bnd_temp(PARA_DATA *para, REAL **var, int var_type, REAL *psi)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN],*af = var[AF], *ab = var[AB];
  REAL *x = var[X],*y = var[Y], *z = var[Z];
  int caseID = para->solv->caseID;

  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;

  /* --------------------------------------------------------------------------
  | wall
  -------------------------------------------------------------------------- */

  j=1;
  FOR_KI
  if(k<=k3)  psi[IX(i,j-1,k)]= para->bc->T_bcS; /*south wall*/
  else psi[IX(i,j-1,k)]= para->bc->T_in;   /*inlet*/
  END_FOR

  j=jmax;
  FOR_KI
  if(k>k1)  psi[IX(i,j+1,k)]= para->bc->T_bcN; /*north wall*/
  else 
  {	  
      psi[IX(i,j+1,k)]= psi[IX(i,j,k)]; /*outflow*/
      an[IX(i,j,k)]=0;
  }
  END_FOR
  
  i=1;
  FOR_JK
  psi[IX(i-1,j,k)]=para->bc->T_bcW;
  END_FOR

  i=imax;
  FOR_JK
  psi[IX(i+1,j,k)]=para->bc->T_bcE;
  END_FOR
 
  k=1;
  FOR_IJ
  psi[IX(i,j,k-1)]=para->bc->T_bcB;
  END_FOR

  k=kmax;
  FOR_IJ
  psi[IX(i,j,k+1)]=para->bc->T_bcT;
  END_FOR

  /* --------------------------------------------------------------------------
  | BOX
  -------------------------------------------------------------------------- */
   j = j1;    
	 for(i=i1+1; i<=i2; i++)
			for(k=1;k<=k2; k++)
	  psi[IX(i,j+1,k)] = para->bc->T_bcBOX;

   j = j2;    
	 for(i=i1+1; i<=i2; i++)
			for(k=1;k<=k2; k++)
	  psi[IX(i,j,k)] = para->bc->T_bcBOX;

  i = i1;  
  for(j=j1+1;j<=j2; j++)
	 for(k=1;k<=k2; k++)
         psi[IX(i+1,j,k)] = para->bc->T_bcBOX;
         
  
  i = i2;  
  for(j=j1+1;j<=j2; j++)
	 for(k=1;k<=k2; k++)
         psi[IX(i,j,k)] = para->bc->T_bcBOX;
 
  k = k2; 
  for(i=i1+1; i<=i2; i++)
	for(j=j1+1;j<=j2; j++)
         psi[IX(i,j,k)] = para->bc->T_bcBOX;

}

/******************************************************************************
|  Set the boundary conditions for pressure
******************************************************************************/
void set_bnd_pressure(PARA_DATA *para, REAL **var, REAL *p)
{
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *af = var[AF], *ab = var[AB];
  int caseID = para->solv->caseID;
  int k1 = para->geom->k1, k2 = para->geom->k2, k3 = para->geom->k3; 
  int i1 = para->geom->i1, i2 = para->geom->i2; 
  int j1 = para->geom->j1, j2 = para->geom->j2;
  
  /* --------------------------------------------------------------------------
  | South Boundary
  ---------------------------------------------------------------------------*/
  j = 1;  
  FOR_KI

	 p[IX(i,j-1,k)] = p[IX(i,j,k)]; as[IX(i,j,k)] = 0.0; 
	 
  END_FOR
  
  /* --------------------------------------------------------------------------
  | North Boundary
  ---------------------------------------------------------------------------*/
  j = jmax;  
   
  FOR_KI
     // if(k<=k1)  p[IX(i,j+1,k)]=0.0;
	 // else {p[IX(i,j+1,k)] = p[IX(i,j,k)]; an[IX(i,j,k)] = 0.0; }
	 p[IX(i,j+1,k)] = p[IX(i,j,k)]; an[IX(i,j,k)] = 0.0;
      //else p[IX(i,j+1,k)]=0.0;
  END_FOR
       
  /* --------------------------------------------------------------------------
  | West Boundary
  -------------------------------------------------------------------------- */
  i = 1;
  FOR_JK
	  p[IX(i-1,j,k)] = p[IX(i,j,k)]; aw[IX(i,j,k)] = 0.0;
  END_FOR
      
  /* --------------------------------------------------------------------------
  | East Boundary
  -------------------------------------------------------------------------- */
  i = imax;      
      FOR_JK
		 p[IX(i+1,j,k)] = p[IX(i,j,k)];  ae[IX(i,j,k)] = 0.0; 
	  END_FOR
   
  // End of East Boundary

    /* --------------------------------------------------------------------------
  | Back Boundary
  ---------------------------------------------------------------------------*/
  k = 1;
     FOR_IJ
	 p[IX(i,j,k-1)] = p[IX(i,j,k)]; ab[IX(i,j,k)] = 0.0;
     END_FOR


  /* --------------------------------------------------------------------------
  | Front Boundary
  ---------------------------------------------------------------------------*/
    k = kmax;  
     FOR_IJ
      p[IX(i,j,k+1)] = p[IX(i,j,k)]; af[IX(i,j,k)] = 0.0;
     END_FOR



  /* --------------------------------------------------------------------------
  | BOX SOUTH
  ---------------------------------------------------------------------------*/

  j=j1;
	 	for(i=i1+1; i<=i2; i++)
		    for(k=1;k<=k2; k++)
			{ p[IX(i,j+1,k)] = p[IX(i,j,k)]; an[IX(i,j,k)] = 0.0;}
  j=j2;
	 	for(i=i1+1; i<=i2; i++)
		    for(k=1;k<=k2; k++)
			{ p[IX(i,j,k)] = p[IX(i,j+1,k)]; as[IX(i,j+1,k)] = 0.0;}

  i=i1;
	 	for(j=j1+1; j<=j2; j++)
		    for(k=1;k<=k2; k++)
			{ p[IX(i+1,j,k)] = p[IX(i,j,k)]; ae[IX(i,j,k)] = 0.0;}

  i=i2;
	 	for(j=j1+1; j<=j2; j++)
		    for(k=1;k<=k2; k++)
			{ p[IX(i,j,k)] = p[IX(i+1,j,k)]; aw[IX(i,j+1,k)] = 0.0;}

  k=k2;
        for(i=i1+1; i<=i2; i++)
	 	   for(j=j1+1; j<=j2; j++)
			{ p[IX(i,j,k)] = p[IX(i,j,k+1)]; ab[IX(i,j,k+1)] = 0.0;}

} // End of set_bnd_pressure( )


/******************************************************************************
| Check the mass conservation of the domain
******************************************************************************/
void mass_conservation(PARA_DATA *para, REAL **var)
{
  int i, j, k;
  int k1= para->geom->k1, k3= para->geom->k3;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *x = var[X], *y = var[Y], *z = var[Z];
  REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL mass_in = 0.0, mass_out = 0.00000001f, mass_ratio;
  REAL sxy,syz,sxz;
  int turn_on = 0;

  /*---------------------------------------------------------------------------
  | Compute the total inflow
  ---------------------------------------------------------------------------*/
  j=0;
   FOR_KI 
	  if(k>k3)   mass_in +=  v[IX(i,j,k)] * (gx[IX(i,j,k)]-gx[IX(i-1,j,k)])*(gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
	  else       mass_in += 0;
   END_FOR
 
  /*---------------------------------------------------------------------------
  | Compute the total outflow
  ---------------------------------------------------------------------------*/

  j=jmax;

   FOR_KI

   if(k<=k1)
   {
   
   sxz= (gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
   sxy= (gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gy[IX(i,j,k)]-gy[IX(i,j-1,k)]);
   syz= (gy[IX(i,j,k)]-gy[IX(i,j-1,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);

   v[IX(i,j,k)]=  v[IX(i,j-1,k)] 
                + (u[IX(i-1,j,k)]- u[IX(i,j ,k)])* syz/sxz
				+ (w[IX(i,j,k-1)]- w[IX(i,j ,k)])* sxy/sxz;
   }

   END_FOR

    j = jmax;
    FOR_KI
	//  if(k<=k1 && v[IX(i,j,k)]>0 )  mass_out +=   v[IX(i,j,k)] * (gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
	//  else       mass_out +=0;

	  	  if(k<=k1  )  mass_out +=   v[IX(i,j,k)] * (gx[IX(i,j,k)]-gx[IX(i-1,j,k)])* (gz[IX(i,j,k)]-gz[IX(i,j,k-1)]);
	  else       mass_out +=0;
    END_FOR
   
  /*---------------------------------------------------------------------------
  | Return the ratio of inflow and outflow
  ---------------------------------------------------------------------------*/
    mass_ratio = mass_in / mass_out;
    printf("mass_in = %f, mass_out = %f, mass_in/mass_out=%f\n", mass_in, mass_out, mass_ratio);

    /*---------------------------------------------------------------------------
    | Adjust the outflow
    ---------------------------------------------------------------------------*/

      j = jmax;
      FOR_KI
		if(k<=k1 )
		{
			//if(v[IX(i,j,k)]>0)  v[IX(i,j,k)] = v[IX(i,j,k)]*mass_ratio;
			//else v[IX(i,j,k)]=0;

			v[IX(i,j,k)] = v[IX(i,j,k)]*mass_ratio;
		}
      END_FOR
   

} // End of mass_conservation()