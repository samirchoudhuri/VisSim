// KEYWORD: EOR2FITS || Converts EOR correlator output to  UV FITS data 
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <ctype.h>
# include <fitsio.h>
# include <stdlib.h>
# include "visfits.h"

char SCANFILE[128], SOURFILE[128];
int NT;

int main(int argc, char *argv[]){

  //timing variables
  int tdiff;
  time_t tstart, tend;

  if(argc!=8){ 
    printf("Usage: %s <input file> <scan file> <source file> <input image  fits file>\n\t <point source file> <output fits file> <No of Threads>\n", argv[0]);
    return 1;
  }
  NT = atoi(argv[7]);
  
  Legals();
  
  CHKFILE(argv[1]);
  
  sscanf(argv[2], "%s", SCANFILE);
  sscanf(argv[3], "%s", SOURFILE);


  CHKFILE(SCANFILE);
  CHKFILE(SOURFILE);
  CHKFILE(argv[4]);
  CHKFILE(argv[5]);

  Read_Inputs(argv[1]); //chkd
  Init_PARM(); 
  
  time (&tstart);
  
  /* Intialize fits file header */
  Init_HDR(argv[6]);
  
  /* Reading iamge fits file */  
  if( Read_Image_FITS(argv[4]) !=0){
    printf("\nERROR: Failed at Read_Image_FITS\n");
    return 1;
  }
  /* Reading point source file */  
  if( Read_Point_Sources(argv[5]) !=0){
    printf("\nERROR: Failed at Read_Point_Source\n");
    return 1;
  }
  /* Write data */
  Write_FITS_data(argv[6]);
  Write_AG_TABLE(argv[6]);
  Write_FQ_TABLE(argv[6]);
  Write_SU_TABLE(argv[6]);
  printf("All extention tables written....\n");

  time(&tend);
  tdiff = (int)difftime (tend, tstart);
  printf("\nYahoooo.... We are DONE\n");
  printf("OUTPUT FILE : %s \n", argv[6]);
  printf("Total runtime  %2d H %2d M %2d S\n", tdiff/3600, tdiff/60, tdiff%60);
  printf("HAVE FUN WITH THE DATA\n");
  printf("=============================================================\n\n\n");

  return 0;
}
