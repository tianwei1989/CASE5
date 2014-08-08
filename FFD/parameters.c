#include "data_structure.h"

/******************************************************************************
|  Input Parameters
******************************************************************************/
void input_para(PARA_DATA *para)
{
  // 1: Lid-Driven Cavity
  // 2: Channel Flow
  // 3: Heat Diffusion in a Cavity
  // 4: Forced Convection
  // 5: Natural Convection in a Square Lid
  // 6: Natural Convection in a tall cavity
  // 7: Mixed Convection in an empty room
  para->solv->caseID = 2;

  /*---------------------------------------------------------------------------
  Initialize the variables
  ---------------------------------------------------------------------------*/
  switch(para->solv->caseID)
  {
    // Miao case 1
    case 1:
      para->geom->imax = 20; //20
      para->geom->jmax = 20; //50
      para->geom->kmax = 20; //20
	  para->geom->Lx   = 2.44f; //3m
      para->geom->Ly   = 2.44f; //m
      para->geom->Lz   = 2.44f; //m
      para->geom->z1   = 0.08f;
	  para->geom->z2   = 2.41f;
	  para->geom->k1   = 4; //10
	  para->geom->k2   = 16;

      para->geom->uniform = 1; //1: uniform; 0: non-uniform 

      para->mytime->dt = 0.1f; 
      para->mytime->t_steady = 100.0f; 
	  para->mytime->t_output =8010;

      para->prob->nu    = 1.53e-5;//1.53e-5f;
      para->prob->tur_model = LAM; //LAM, CHEN, CONST(==101nu)
      para->prob->chen_a = 0.03874f;  /* the coeffcient of Chen's model*/
      para->solv->solver = GS;
      para->solv->advection_solver = SEMI;  //LAX, SEMI, UPWIND, UPWIND_NEW;
      para->solv->interpolation = BILINEAR; //BILINEAR, FSJ
      para->solv->check_residual = 0; 
      para->solv->read_file = 0; 

      para->outp->version = DEBUG; //DEMO, DEBUG;

      para->bc->bcN = OUTFLOW;
      para->bc->bcS = INFLOW;
      para->bc->bcW = NOSLIP;
      para->bc->bcE = NOSLIP;   
      para->bc->bcF = NOSLIP;
      para->bc->bcB = NOSLIP;

      para->bc->VY_bcS = 1.36f;

      break;   
      
    // Miao 3
    case 2:
      para->geom->imax = 20;
      para->geom->jmax = 20;
      para->geom->kmax = 20;
      para->geom->Lx   = 2.44;
      para->geom->Ly   = 2.44;
      para->geom->Lz   = 2.44;
	  para->geom->z1   = 0.08f;
	  para->geom->z2   = 1.22f;
	  para->geom->z3   = 2.41f;
	  para->geom->k1   = 3; //10
	  para->geom->k2   = 10;
	  para->geom->k3   = 17;
	  para->geom->x1   = 0.61;
	  para->geom->x2   = 1.83;
      para->geom->i1   = 6; //10
	  para->geom->i2   = 14;
	  para->geom->y1   = 0.61;
	  para->geom->y2   = 1.83;
      para->geom->j1   = 6; //10
	  para->geom->j2   = 14;
      para->geom->uniform = 1; //1: uniform; 0: non-uniform 

      para->mytime->dt = 0.1; 
      para->mytime->t_output = 1000;
      para->mytime->t_steady = 490; 
      para->prob->nu = 1.789e-5;//1.53e-5f; 1.53e-5;
      para->prob->tratio=1.0;
      para->prob->tur_model = LAM; //LAM, CHEN, CONST(==101nu)
      para->prob->chen_a = 0.03874;  /* the coeffcient of Chen's model*/
	  para->prob->alpha =1.96e-5;
	  para->prob->beta =0.0;//0.0032; //0.0032;
	  para->prob->gravx =0;
	  para->prob->gravy =0;
	  para->prob->gravz =-9.8;
	  para->prob->Temp_opt= 288.16;
   
	  para->outp->version = DEBUG; //DEMO, DEBUG;

      para->solv->solver = GS; // GS, TDMA
      para->solv->advection_solver = SEMI;  //LAX, SEMI, UPWIND, UPWIND_NEW;
      para->solv->interpolation = BILINEAR; //BILINEAR, FSJ
      para->solv->check_residual = 0;

      para->bc->bcN = OUTFLOW;        
      para->bc->bcS = INFLOW;
      para->bc->bcW = NOSLIP;
      para->bc->bcE = NOSLIP; 
      para->bc->bcF = NOSLIP;
      para->bc->bcB = NOSLIP;
      
      para->bc->VY_bcS = 1.36;

      para->bc->T_bcN = 300.86f;
      para->bc->T_bcS = 300.77f;
      para->bc->T_bcW = 300.1f;
      para->bc->T_bcE = 300.35f;
	  para->bc->T_bcT = 298.98f;
      para->bc->T_bcB = 300.13f;
      para->bc->T_in  = 295.15f;
      para->bc->T_bcBOX = 308.3f;

      break;   

    // Heat Diffusion in a Cavity
    case 3:
      para->geom->imax = 15;
      para->geom->jmax = 15;
      para->geom->Lx   = 1.0;
      para->geom->Ly   = 1.0;             
      para->geom->uniform = 1; //1: uniform; 0: non-uniform 

      para->mytime->dt = 1.0f; 
      para->mytime->t_output = 100;
      para->mytime->t_steady = 900.0f; 
      para->prob->alpha = 0.001f; 
      para->prob->mu    = 1.0f;
      para->prob->rho   = 10000.0f;        

      para->bc->bcN = NOSLIP;
      para->bc->bcS = NOSLIP;
      para->bc->bcW = INFLOW;
      para->bc->bcE = OUTFLOW;   

      para->bc->VX_bcN = 0.0f;
      para->bc->VX_bcS = 0.0f;
      para->bc->VX_bcE = 0.0f;
      para->bc->VX_bcW = 0.0f;

      para->bc->VY_bcN = 0.0f;
      para->bc->VY_bcS = 0.0f;
      para->bc->VY_bcE = 0.0f;
      para->bc->VY_bcW = 0.0f;

      para->bc->bcTW = TCONST;
      para->bc->bcTE = TCONST;
      para->bc->bcTS = TCONST;
      para->bc->bcTN = TCONST;   

      para->bc->T_bcN = 1.0f;
      para->bc->T_bcS = 1.0f;
      para->bc->T_bcW = 0.0f;
      para->bc->T_bcE = 0.0f;

      break;  
    /*-------------------------------------------------------------------------
    | Forced Convection
    -------------------------------------------------------------------------*/
    case 4:
      para->geom->imax = 143;
      para->geom->jmax = 71;
      para->geom->j1 = 32;
      para->geom->j2 = 56;
      para->geom->y1 = 0.48;
      para->geom->y2 = 0.168;

      para->geom->Lx   = 9.0f;
      para->geom->Ly   = 3.0f;             
      para->geom->uniform = 0; //1: uniform; 0: non-uniform 

      para->mytime->dt = 0.005f; 
      para->mytime->t_output = 220000;
      para->mytime->t_steady = 1000.0f; 
      para->prob->mu = 15.3e-6f; 
      para->prob->rho = 1.0f;

      para->prob->tur_model = LAM; //LAM, CHEN, CONST(==101nu)
      para->prob->chen_a = 0.03874f;  /* the coeffcient of Chen's model*/

      para->solv->solver = GS;
      para->solv->advection_solver = UPWIND;  //LAX, SEMI, UPWIND, UPWIND_NEW;
      para->solv->interpolation = BILINEAR; //BILINEAR, FSJ
      para->solv->check_residual = 0;

      para->bc->bcN = NOSLIP;
      para->bc->bcS = NOSLIP;
      para->bc->bcW = INFLOW;
      para->bc->bcE = OUTFLOW;   
      
      para->outp->version = DEMO; //DEMO, DEBUG;

      para->bc->VX_bcW = 0.455f;
      break;  
    /*-------------------------------------------------------------------------
    | Natural Convection in a Square Lid
    -------------------------------------------------------------------------*/
    case 5:
      para->geom->imax = 15;
      para->geom->jmax = 15;
      para->geom->Lx   = 1.0f;
      para->geom->Ly   = 1.0f;             
      para->geom->uniform = 1; //1: uniform; 0: non-uniform 

    	para->mytime->dt = 0.1f; 
      para->mytime->t_output = 1000;
      para->mytime->t_steady = 8000.0f; 

      para->prob->nu = 1.0f/2800.0f;
      para->prob->alpha = 0.001f;      
      para->prob->tur_model = LAM; //LAM, CHEN, CONST(==101nu)
      para->prob->chen_a = 0.03874f;  /* the coeffcient of Chen's model*/
      para->prob->gravz = -9.8f;
      para->prob->beta = 3.186e-3f;

      para->solv->solver = GS; // GS, TDMA
      para->solv->advection_solver = SEMI;  //LAX, SEMI, UPWIND, UPWIND_NEW;
      para->solv->interpolation = BILINEAR; //BILINEAR, FSJ
      para->solv->check_residual = 0;

      para->bc->bcN = NOSLIP;        
      para->bc->bcS = NOSLIP;
      para->bc->bcW = NOSLIP;
      para->bc->bcE = NOSLIP;   
      
      para->bc->bcTW = TCONST;
      para->bc->bcTE = TCONST;
      para->bc->bcTS = ADIBATIC;
      para->bc->bcTN = ADIBATIC;   

      para->bc->T_bcW = 0.0f;
      para->bc->T_bcE = 1.0f;      
      break;   
    /*-------------------------------------------------------------------------
    | Natural Convection in a Tall Cavity
    -------------------------------------------------------------------------*/
    case 6:
      para->geom->imax = 19;
      para->geom->jmax = 39;
      para->geom->Lx   = 0.076f;
      para->geom->Ly   = 2.18f;             
      para->geom->uniform = 0; //1: uniform; 0: non-uniform 

    	para->mytime->dt = 0.01f; 
      para->mytime->t_output = 50000;
      para->mytime->t_steady = 400.0f; 

      para->prob->mu = 1.6956e-5f;
      para->prob->rho = 1.1614f;
      para->prob->alpha = 2.074e-5f;      
      para->prob->tur_model = LAM; //LAM, CHEN, CONST(==101nu)
      para->prob->chen_a = 0.03874f;  /* the coeffcient of Chen's model*/
      para->prob->gravz = -9.8f;
      para->prob->beta = 3.47e-3f;

      para->solv->solver = GS; // GS, TDMA
      para->solv->advection_solver = SEMI;  //LAX, SEMI, UPWIND, UPWIND_NEW;
      para->solv->interpolation = BILINEAR; //BILINEAR, FSJ
      para->solv->check_residual = 0;

      para->bc->bcN = NOSLIP;        
      para->bc->bcS = NOSLIP;
      para->bc->bcW = NOSLIP;
      para->bc->bcE = NOSLIP;   
      
      para->bc->bcTW = TCONST;
      para->bc->bcTE = TCONST;
      para->bc->bcTS = ADIBATIC;
      para->bc->bcTN = ADIBATIC;   

      para->bc->T_bcW = 0.0f;
      para->bc->T_bcE = 19.6f;      
      break;   
    /*-------------------------------------------------------------------------
    | Mixed Convection in an empty room
    -------------------------------------------------------------------------*/
    case 7:
      para->geom->imax = 30;
      para->geom->jmax = 30;
      para->geom->Lx   = 1.04f;
      para->geom->Ly   = 1.04f;   
      para->geom->j1 = 6;
      para->geom->j2 = 24;
      para->geom->y1 = 0.024;
      para->geom->y2 = 0.018;
      para->geom->uniform = 0; //1: uniform; 0: non-uniform 

    	para->mytime->dt = 0.1f; 
      para->mytime->t_output = 5000;
      para->mytime->t_steady = 490.0f; 

      para->outp->version = DEBUG; //DEMO, DEBUG;

      // Air at 298.4K=25.25C
      para->prob->mu = 1.812e-5f;
      para->prob->rho = 1.1861f;
      para->prob->alpha = 2.376846e-5f;      
      para->prob->tur_model = LAM; //LAM, CHEN, CONST(==101nu)
      para->prob->chen_a = 0.03874f;  /* the coeffcient of Chen's model*/
      para->prob->gravz = -9.8f;
      para->prob->beta = 3.43e-3f;

      para->solv->solver = GS; // GS, TDMA
      para->solv->advection_solver = SEMI;  //LAX, SEMI, UPWIND, UPWIND_NEW;
      para->solv->interpolation = BILINEAR; //BILINEAR, FSJ
      para->solv->check_residual = 0;

      para->bc->bcN = NOSLIP;        
      para->bc->bcS = NOSLIP;
      para->bc->bcW = INFLOW;
      para->bc->bcE = OUTFLOW;   
      para->bc->VX_bcW = 0.57f;
      
      para->bc->bcTW = TCONST;
      para->bc->bcTE = TCONST;
      para->bc->bcTS = TCONST;
      para->bc->bcTN = TCONST;   
      
      para->bc->T_bcW = 15.0f;
      para->bc->T_bcE = 15.0f;      
      para->bc->T_bcN = 15.0f;      
      para->bc->T_bcS = 35.5f;      

      para->outp->Temp_ref = 25.25;
      para->prob->diff = 0.00001;
      para->prob->force = 0.1; 
      para->prob->source = 0.1;
      para->outp->winx   = 600;
      para->outp->winy   = 600;
      para->outp->v_ref  = 1.0; 
      break;   
  }


} // End of input_para( )