# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <ctype.h>
# include <fitsio.h>
# include <stdlib.h>
# include <unistd.h>
# include "eor2fits.h"


void CHKFILE(char filename[]){

  if(access(filename, F_OK)!=0){
    printf("File %s does not exists\n", filename);
    printf("Exiting to system\n\n");
    exit(0);
  }
}
void COPYLEFT(){

  char param[128];
  FILE *fp;
  CHKFILE(COPYLEFTFILE);
  fp = fopen(COPYLEFTFILE, "r");
  while(feof(fp)==0){
    if(fgets(param,100,fp),(param[0]=='$'))
      printf("%s", param);
  }
  fclose(fp);
}

void replace_nulls(char *str, int len){ 
int i,k ;
 k = strlen(str) ;
 for (i=k; i<len; i++) str[i] = ' ' ;
}
double jd2gst(double jd){
  // gst in degrees corresponding to jd
  // adopted from gvfits check this ****
  double c0,c1,c2,c3;
  double temp,t;
  c0=24110.54841;
  c1=8640184.812866;
  c2=0.093104;
  c3=-6.2E-6;
  t=(jd-2451545.0)/36525;
  temp=((c0+c1*t+c2*t*t+c3*t*t*t))/240.0;
  return (temp);
}
char *jd2iau_date(double J){
  
  // This function is adopted from GMRT gvfits file "csubs1.c" and modified 
  // to meet necessity here.
 
  static char date[32] ;
  int month, day;
  long year, a, c, d, x, y, jd;
  double dd;
  
  if( J < 1721425.5 ) return 0 ; /* January 1.0, 1 A.D.  dont accept BC ! */
  
  jd = J + 0.5; /* round Julian date up to integer */
  
  /* Find the number of Gregorian centuries
   * since March 1, 4801 B.C.
   */
  a = (100*jd + 3204500L)/3652425L;
  
  /* Transform to Julian calendar by adding in Gregorian century years
   * that are not leap years.
   * Subtract 97 days to shift origin of JD to March 1.
   * Add 122 days for magic arithmetic algorithm.
   * Add four years to ensure the first leap year is detected.
   */
  c = jd + 1486;
  if( jd >= 2299160.5 )
    c += a - a/4;
  else
    c += 38;
  /* Offset 122 days, which is where the magic arithmetic
   * month formula sequence starts (March 1 = 4 * 30.6 = 122.4).
   */
  d = (100*c - 12210L)/36525L;
  x = (36525L * d)/100L; /* Days in that many whole Julian years */
  
  /* Find month and day. */
  y = ((c-x)*100L)/3061L;
  day = c - x - ((306L*y)/10L);
  month = y - 1;
  if( y > 13 ) month -= 12;
  
  /* Get the year right. */
  year = d - 4715;
  if( month > 2 ) year -= 1;
  
  a = (jd + 1) % 7; /* Day of the week. */
  
  //    dd = day + J - jd + 0.5; /* Fractional part of day. */
  //    day = dd;
  
  sprintf(date,"%04ld-%02d-%02d / UT",year,month,day) ;
  
  return date ;
}


double HMS2D(float HH, float MM, float SS){

  return (HH + MM/60. + SS/3600.)*15.;
}

