#include <stdio.h>
#include <string.h>

#include "data_structure.h"
#include "read_data.h"

FILE *file_params;

/******************************************************************************
| Write the data to a file for Tecplot 
******************************************************************************/
int read_data(PARA_DATA *para, REAL **var)
{
  int i,j, k;
  int tempi, tempj, tempk;
  int imax = para->geom->imax;
  int jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  char *temp, string[400];

  if( (file_params=fopen("steady.plt","r")) == NULL ) 			
  {
    fprintf(stderr,"Error:can not open error file!\n");
    return 0;
  }
    
  temp = fgets(string, 400, file_params);
  temp = fgets(string, 400, file_params);
  temp = fgets(string, 400, file_params);    

  FOR_ALL_CELL
    temp = fgets(string, 400, file_params);
    sscanf( string,"%f%f%f%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
       &var[X][IX(i,j,k)], &var[Y][IX(i,j,k)], &var[Z][IX(i,j,k)], &tempi, &tempj, &tempk,  
       &var[VXM][IX(i,j,k)], &var[VYM][IX(i,j,k)], &var[VZM][IX(i,j,k)], &var[TEMPM][IX(i,j,k)], 
       &var[VX][IX(i,j,k)], &var[VY][IX(i,j,k)], &var[VZ][IX(i,j,k)], &var[TEMP][IX(i,j,k)],
       &var[DEN][IX(i,j,k)], &var[IP][IX(i,j,k)], 
       &var[AP][IX(i,j,k)], &var[AW][IX(i,j,k)], &var[AE][IX(i,j,k)], 
       &var[AN][IX(i,j,k)], &var[AS][IX(i,j,k)], &var[AF][IX(i,j,k)],
       &var[AB][IX(i,j,k)], &var[B][IX(i,j,k)], 
       &var[TMP1][IX(i,j,k)], &var[TMP2][IX(i,j,k)], &var[TMP3][IX(i,j,k)]);
  END_FOR      
    
  fclose(file_params);

  return 1;

} // End of read_dara()

