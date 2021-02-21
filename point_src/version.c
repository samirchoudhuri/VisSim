# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include "version.h"
# include "eor2fits.h"
char versionstr[7];
char get_version(){
  char hversion[2], lversion[2];
  if((int)HVERSION==0)
    sprintf(hversion, "00");
  else  if((int)HVERSION<10)
    sprintf(hversion, "0%1d", (int)HVERSION);
  else
    sprintf(hversion, "%2d", (int)HVERSION);
  
  if(LVERSION==0)
    sprintf(lversion, "00");
  else  if((int)LVERSION<10)
    sprintf(lversion, "0%1d", (int)LVERSION);
  else
    sprintf(lversion, "%2d", (int)LVERSION);
  
  sprintf(versionstr,"%2s%2d.%2s", hversion, (int)IATUTC, lversion);
}

int main(int argc, char *argv[]){
  
  get_version();
  printf("%s\n", versionstr);
  return 0;
}
