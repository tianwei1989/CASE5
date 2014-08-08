#include <stdio.h>
#include "data_structure.h"
#include "solver.h"
#include "write_data.h"
#include "diffusion.h"
#include "projection.h"
#include "advection.h"
#include "timing.h"
#include "solver_gs.h"
#include "solver_tdma.h"
#include "boundary.h"
#include "utility.h"

/******************************************************************************
| FFD Solver
******************************************************************************/
void FFD_solver(PARA_DATA *para, REAL **var)
{
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int size = (imax+2) * (jmax+2) * (kmax+2);
  int t_step = 0, t_output = para->mytime->t_output;
  REAL t_steady = para->mytime->t_steady;
  REAL dt = para->mytime->dt;
  REAL *u = var[VX], *v = var[VY], *w = var[VZ];
  REAL *den = var[DEN], *temp = var[TEMP];
  REAL *u_mean = var[VXM], *v_mean = var[VYM], *w_mean = var[VZM];
  REAL *temp_mean = var[TEMPM];
 
  int cal_mean = para->outp->cal_mean;
  int i;

  /*---------------------------------------------------------------------------
  | Solver Loop
  ---------------------------------------------------------------------------*/
  while( para->mytime->t_step < t_output)
  {
    vel_step(para, var);	  
    temp_step(para, var);
    //den_step(para, var);

    timing(para);

    if(para->mytime->t>t_steady && cal_mean==0)
    {
      cal_mean = 1;
	  t_step += 1;
      printf("start to calculate mean properties.\n");
    }   

   if(cal_mean == 1)
     for(i=0; i<size; i++)
     {
       u_mean[i] += u[i];
       v_mean[i] += v[i];
       w_mean[i] += w[i];
       temp_mean[i] += temp[i];
	   
     }
  } // End of While loop  

  /*---------------------------------------------------------------------------
  | Post Process
  ---------------------------------------------------------------------------*/
  if(cal_mean == 1)
    for(i=0; i<=size; i++)
    {
      u_mean[i] = u_mean[i] / t_step;
      v_mean[i] = v_mean[i] / t_step;
      w_mean[i] = w_mean[i] / t_step;
      temp_mean[i] = temp_mean[i] / t_step;
    }

  write_data(para, var, "result");
  //write_data1(para, var, "result1");
  para->prob->output = 1;
} // End of FFD_solver( ) 


/******************************************************************************
| calculate the temperature
******************************************************************************/ 
void temp_step(PARA_DATA *para, REAL **var)
{
  REAL *T = var[TEMP], *T0 = var[TMP1];  	  
  
  advection(para, var, TEMP, T0, T); 	   	
  diffusion(para, var, TEMP, T, T0);

} // End of temp_step( )

/******************************************************************************
| calculate the temperature
******************************************************************************/ 
void den_step(PARA_DATA *para, REAL **var)
{
  REAL *den = var[DEN], *den0 = var[TMP1];
  	  
  advection(para, var, DEN, den0, den); 	   	
  diffusion(para, var, DEN, den, den0); 	

} // End of den_step( )

/******************************************************************************
| calculate the velocity
******************************************************************************/ 
void vel_step(PARA_DATA *para, REAL **var)
{
  REAL *u  = var[VX],  *v  = var[VY],    *w  = var[VZ];
  REAL *u0 = var[TMP1], *v0 = var[TMP2], *w0 = var[TMP3];
  	  

  advection(para, var, VX, u0, u);
  advection(para, var, VY, v0, v); 	    
  advection(para, var, VZ, w0, w); 	    

 diffusion(para, var, VX, u, u0);   
  diffusion(para, var, VY, v, v0); 
 diffusion(para, var, VZ, w, w0); 

 mass_conservation(para, var);
  projection(para, var); 


  //swap(para, var);

  //mass_conservation(para, var);

  //printf("finish the advection \n");
 

  //projection(para, var);    
  //mass_conservation(para, var);
    
} // End of vel_step( )


/******************************************************************************
| Solver of equations
******************************************************************************/ 
void equ_solver(PARA_DATA *para, REAL **var, int var_type, REAL *psi)
{
  switch(para->solv->solver)
  {
    case GS:
      Gauss_Seidel(para, var, var_type, psi);
      break;
    case TDMA:
      TDMA_3D(para, var, var_type, psi);
      break;
  }
}// end of equ_solver