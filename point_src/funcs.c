// KEYWORD: UVFITS_FUNCS || functions to simulate UV FITS data 
// COMMENTS: 1) chk what are the parameters written by write group hdr
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <ctype.h>
# include <fitsio.h>
# include <stdlib.h>
# include <unistd.h>
# include <omp.h>
# include <nr.h>//only for julday
# include <nrutil.h>
# include <time.h>
# include "visfits.h"
# include "version.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
double Beam(double theta,double freq);
// INPUTS
char   TELESC[128], ANTFILE[128];
float  TELE_LAT_DD, TELE_LAT_Mm, TELE_LAT_SS;
float  TELE_LON_DD, TELE_LON_Mm, TELE_LON_SS;

float  FREQOBS, BWOBS;
double RESTFREQ;

char   NAMEOBS[8];
int    NCHAN;
float  INTIME;
int    OBS_DD, OBS_MM, OBS_YY;
float  ST_HH, ST_MM, ST_SS;

// for random number generation
unsigned long int seedA, seedP, seedN,seedPos,seedran;
gsl_rng *rA,*rP,*rN,*rPos,*rran;
double  SIG_A, SIG_P, SIG_N, SIG_Acont, SIG_Pcont,SIG_Pos;
// done 

int    NTbin,tmpbin;
double xshift=0.,yshift=0.;
double Uval,Umax,flagran; 
// INCODE parameters
int    NSOUR;
long   NSCAN, NBASELINE, group, GCOUNT, NRECORD;
float  CWOBS;

float  TELE_LAT, TELE_LON;
double DATEOBS, STTIM;
char   DATEOBS_S[16];
float *souRA, *souDC;

double  *bx, *by, *bz;
double  *u, *v, *w;
double  LTT;

long VISlen, npixels; 
float *VISword, *randpar, *Idata;
float *xPS, *yPS,*xPS0, *yPS0, **IPS,*alphapt;
float *gA_re, *gA_im;

int    NGridI, NPSour;
float  DIx, DIy;

extern char SCANFILE[128], SOURFILE[128];
extern int NT;
char versionstr[10];

void get_version(){
    sprintf(versionstr,"%02d%2d.%02d", HVERSION, (int)IATUTC, LVERSION);
}

void printerror(int status){
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
  
  if (status){
      fits_report_error(stderr, status); // Print error report
      exit( status );    // Terminate the program, returning error status
    }
}
void Calc_uvw(double hh, int SUID){
 
  int ii;
  float dec;
  // hh is input in degrees
  // dec is input in degrees

  hh = hh*DR;
  dec = souDC[SUID-1]*DR;

  for(ii=0;ii<NANTE;ii++)
    {
      u[ii] = sin(hh)*bx[ii] + cos(hh)*by[ii] ;
      v[ii] = -sin(dec)*cos(hh)*bx[ii] + sin(dec)*sin(hh)*by[ii] + cos(dec)*bz[ii];
      w[ii] =  cos(dec)*cos(hh)*bx[ii] - cos(dec)*sin(hh)*by[ii] + sin(dec)*bz[ii];
    }
}

void Legals(){

  time_t runtime;

  COPYLEFT();
  get_version();
  
  time(&runtime);
  printf("\n$ ************************************************************* $\n");
  printf("  \t\t RUNNING eor2fits VERSION: %7s \n", versionstr);
  printf("  \t\t TIME: %s", ctime(&runtime));
  printf("$ ************************************************************* $\n\n");
  
}
void Read_Inputs(char *inpfile){

  FILE *fp;
  char param[128];

  CHKFILE(inpfile);

  fp=fopen(inpfile,"r");

  // TELESCOPE PARAMETERS (4)
  
  fgets(param,100, fp);
  fgets(param,100, fp);
  fgets(param,100, fp);

  fgets(param,100, fp);
  sscanf(param,"%s", TELESC);
  
  fgets(param,100, fp);
  sscanf(param,"%s", ANTFILE);
  
  fgets(param,100, fp);
  sscanf(param,"%e%e%e", &TELE_LAT_DD, &TELE_LAT_Mm, &TELE_LAT_SS);
  
  fgets(param,100, fp);
  sscanf(param,"%e%e%e", &TELE_LON_DD, &TELE_LON_Mm, &TELE_LON_SS);
  
  // FREQUENCY PARAMETERS (3)
  fgets(param,100, fp);
  fgets(param,100, fp);
  fgets(param,100, fp);
  fgets(param,100, fp);

  fgets(param,100, fp);
  sscanf(param,"%e", &FREQOBS);
  
  fgets(param,100, fp);
  sscanf(param,"%e", &BWOBS);

   fgets(param,100, fp);
  sscanf(param,"%lf", &RESTFREQ);

  // OBSERVATION PARAMETERS (4)
  fgets(param,100, fp);
  fgets(param,100, fp);
  fgets(param,100, fp);
  fgets(param,100, fp);

  fgets(param,100, fp);
  sscanf(param,"%s", NAMEOBS);

  fgets(param,100, fp);
  sscanf(param,"%d", &NCHAN);

  fgets(param,100, fp);
  sscanf(param,"%f", &INTIME);

  fgets(param,100, fp);
  sscanf(param,"%d%d%d", &OBS_DD, &OBS_MM, &OBS_YY);

  fgets(param,100, fp);
  sscanf(param,"%f%f%f", &ST_HH, &ST_MM, &ST_SS);


  // for random number generation

  fgets(param,100, fp);
  sscanf(param,"%lf%ld", &SIG_A,&seedA);
  fgets(param,100, fp);
  sscanf(param,"%lf%ld", &SIG_P,&seedP);
  fgets(param,100, fp);
  sscanf(param,"%lf%ld", &SIG_N,&seedN);
  fgets(param,100, fp);
  sscanf(param,"%lf", &SIG_Acont);
  fgets(param,100, fp);
  sscanf(param,"%lf", &SIG_Pcont);
  fgets(param,100, fp);
  sscanf(param,"%d", &NTbin);
  fgets(param,100, fp);
  sscanf(param,"%lf%ld", &SIG_Pos,&seedPos);
  fgets(param,100, fp);
  sscanf(param,"%lf%ld%lf", &Umax,&seedran,&flagran);
  fclose(fp);
  printf("Inputs are read....\n");

  // for random number generation 
  rA= gsl_rng_alloc(gsl_rng_cmrg);
  rP= gsl_rng_alloc(gsl_rng_cmrg);
  rN= gsl_rng_alloc(gsl_rng_cmrg);
  rPos= gsl_rng_alloc(gsl_rng_cmrg);
  rran= gsl_rng_alloc(gsl_rng_cmrg);
  
  gsl_rng_set (rA, seedA);
  gsl_rng_set (rP, seedP);
  gsl_rng_set (rN, seedN);
  gsl_rng_set (rPos, seedPos);
  gsl_rng_set (rran, seedran);

  // done

}
void Init_PARM(){

  int ii, tempi, ERECORD;
  FILE *fp;
  double TJD, JD0;
  char TIMEFILE[128], param[128];

  // print the version info, compilation date and copyright info
  
  printf("Initializing parameters....\n");

  // checking if the files exists
  CHKFILE(ANTFILE);

  fp = fopen(SCANFILE, "r");
  fgets(param,100, fp);
  fgets(param,100, fp);
  
  NSCAN=-1;
  while(feof(fp)==0){
    fscanf(fp, "%*d%*d%*d%d%*d", &tempi);
    NSCAN++;
  }
  fclose(fp);

  printf("Read existance of %ld scans\n", NSCAN);
  if(NSCAN==0){
    printf("ERROR: No Scan info in the scan file\n");
    exit(1)  ;
  }
  int BR, ER;
  NRECORD = 0;
  fp = fopen(SCANFILE, "r");
  fgets(param,100, fp);
  fgets(param,100, fp);
  for(ii=0; ii<NSCAN; ii++){
    fscanf(fp, "%*d%*d%*d%d%d", &BR, &ER);
    NRECORD += (ER-BR+1);
  }
  fclose(fp);

  printf("NRECORD = %ld\n", NRECORD);

  NBASELINE = (long)(NANTE-1)*NANTE/2;
  CWOBS  = (float)BWOBS/(1.*NCHAN);
  GCOUNT = NBASELINE*NRECORD; 

  // setting the julday of observation and TIME1

  DATEOBS = 1.*julday(OBS_MM, OBS_DD, OBS_YY);
  STTIM   = HMS2D(ST_HH, ST_MM, ST_SS);
 
  TELE_LAT = HMS2D(TELE_LAT_DD, TELE_LAT_Mm, TELE_LAT_SS)/15.;
  TELE_LON = HMS2D(TELE_LON_DD, TELE_LON_Mm, TELE_LON_SS)/15.;
 
  // setting Light Travel TIme
  LTT = FREQOBS*1.e6/(1.*C);

  // setting arrays
  bx = (double*)calloc(NANTE, sizeof(double));
  by = (double*)calloc(NANTE, sizeof(double));
  bz = (double*)calloc(NANTE, sizeof(double));

  u = (double*)calloc(NANTE, sizeof(double));
  v = (double*)calloc(NANTE, sizeof(double));
  w = (double*)calloc(NANTE, sizeof(double));
  
  fp = fopen(ANTFILE, "r");
  for(ii=0; ii<NANTE; ii++)
    fscanf(fp, "%*s%lf%lf%lf\n", &bx[ii], &by[ii], &bz[ii]);
  fclose(fp);

  // Reding source info
  float  SOURA_HH, SOURA_Mm, SOURA_SS, SOURA;
  float  SOUDC_DD, SOUDC_Mm, SOUDC_SS, SOUDC;

  fp = fopen(SOURFILE, "r");
  fgets(param,100, fp);
  fgets(param,100, fp);
  
  NSOUR = -1;
  while(feof(fp)==0){
    fscanf(fp, "%*s%d%*s%*f%*f%*f%*f%*f%*f", &tempi);
    NSOUR++;
  }
  fclose(fp);  
  printf("Read existance of %d sources\n", NSOUR);

  if(NSOUR==0){
    printf("ERROR: NO source info in source file\n");
    exit(1);
  }
  //  printf("NSOUR = %d\n", NSOUR);

  fp = fopen(SOURFILE, "r");
  fgets(param,100, fp);
  fgets(param,100, fp);

  souRA = (float*)calloc(NSOUR, sizeof(float));
  souDC = (float*)calloc(NSOUR, sizeof(float));

  for(ii=0; ii<NSOUR; ii++){
    fscanf(fp, "%*s%*d%*s%f%f%f%f%f%f", &SOURA_HH, &SOURA_Mm, &SOURA_SS, 
	   &SOUDC_DD, &SOUDC_Mm, &SOUDC_SS);

    SOURA = HMS2D(SOURA_HH, SOURA_Mm, SOURA_SS);
    SOUDC = HMS2D(SOUDC_DD, SOUDC_Mm, SOUDC_SS)/15.;
   
    souRA[ii] = SOURA;
    souDC[ii] = SOUDC;
  }
  fclose(fp);


  float ST_ha, ED_ha, timestamp, gst;
  timestamp = STTIM;
  gst = jd2gst(DATEOBS + timestamp );

  gst = gst-floor(gst) + 1.*((int)gst)/360.;
  ST_ha  = gst + TELE_LON - souRA[0];
  ED_ha  = ST_ha + 360.*(1.*ER*INTIME)/86400.;
  ST_ha  /= 15.;
  ED_ha  /= 15.;
  printf("gst = %e\n", gst);
  printf("DATEOBS = %f\n", DATEOBS);
  printf("Starting Hour Angle : %.2f\n", ST_ha);
  printf("Finishing Hour Angle : %.2f\n", ED_ha);

  printf("\nParameter values intialized...\n");
}
void Init_HDR(char *FITSFILE){

  //cfitsio INPUTS
  fitsfile *fptr;
  int status;
  int extend, blocked, groups;

  //control inputs
  int ii;
  char OPT[1];

  //deader variables
  long naxes[NAXIS];
  char keynam[6];

  printf("Initialising Header....\n");

  if(access(FITSFILE, F_OK)==0){
    
    printf("\nFile %s already exists\n", FITSFILE);
    printf("Will Exit to system\n\n");
    exit(1);
  }
  // remove(FITSFILE);
  status = 0;
  if(fits_create_file(&fptr, FITSFILE, &status))
    printerror(status);
 
  naxes[0] = 0;	
  naxes[1] = NCMPLX;   
  naxes[2] = NSTOKES;  
  naxes[3] = NCHAN;
  naxes[4] = NSIDE;	
  naxes[5] = 1;
  naxes[6] = 1;
  
  extend   = TRUE;
  blocked  = TRUE;
  groups   = TRUE; 
  
  printf("Now will write header....\n");

  sprintf(DATEOBS_S, "%s", jd2iau_date(DATEOBS));
  
  if(fits_write_grphdr(fptr, TRUE, FLOAT_IMG, NAXIS, naxes, 
		       PCOUNT, GCOUNT, extend, &status)) 
    printerror(status);
  if(fits_write_key_log(fptr, "BLOCKED", TRUE, "Tape may be blocked ", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "OBJECT", "MULTI", "", &status)) 
    printerror(status);
  if(fits_write_key_str(fptr, "TELESCOP", TELESC, "Telescope Name", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "INSTRUME", "SIMUL", "Correlator Name", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "OBSERVER", NAMEOBS, "Observer Name", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "PROJECT", "SIMUL", "Project Name", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "DATE-OBS", DATEOBS_S, "Obs start date YYYY-MM-DD ", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "DATE-MAP", DATEOBS_S, "Last processing date YYYY-MM-DD", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "SORT", "TB", "UV data in time-baseline order", &status))
    printerror(status);
  if(fits_write_key_dbl(fptr, "BSCALE", 1.0E+00, 11, "REAL = TAPE * BSCALE + BZERO", &status))
    printerror(status);
  if(fits_write_key_dbl(fptr, "BZERO", 0.0E+00, 11, " ", &status))
    printerror(status);
  if(fits_write_key_str(fptr, "BUNIT", "UNCALIB", "Units of flux", &status))
    printerror(status);
  if(fits_write_key_dbl(fptr, "EPOCH", EPOACH, 9, "Epoch of RA DEC", &status))
    printerror(status);
  /* if(fits_write_key_lng(fptr, "VELREF", 259, ">256 RADIO, 1 LSR 2 HEL 3 OBS", &status)) */
  /*   printerror(status); */
  // hard coded, find out what
  /* if(fits_write_key_dbl(fptr, "ALTRVAL", 2.3E+05, 11,  "Altenate FREQ/VEL ref pixel", &status)) */
  /*   printerror(status); */
  if(fits_write_key_dbl(fptr, "ALTRPIX", 6.4E+01, 11,  "Altenate FREQ/VEL ref pixel", &status))
    printerror(status);
  /* if(fits_write_key_dbl(fptr, "RESTFREQ", RESTFREQ*1.e6,  11, "Rest FREQ of line", &status)) */
  /*   printerror(status); */

  //declare ctypes
  char *CTYPE[6] =  {"COMPLEX", "STOKES", "FREQ", "IF", "RA", "DEC"};
  /* double  CRVAL[6] =  {1., -1., (FREQOBS-BWOBS/2.+CWOBS/2.)*1.E6, 1., 0., 0.}; */
  /* double  CDELT[6] =  {1., -1., CWOBS*1.E6,   1., 1.,    1.    }; */
  /* double  CRPIX[6] =  {1.,  1., 0.5,     1., 1.,    1.    }; */
  /* double  CROTA[6] =  {0.,  0., 0.,      0., 0.,    0.,   }; */
  double  CRVAL[6] =  {1., -1., (FREQOBS)*1.E6, 1., 0., 0.};
  double  CDELT[6] =  {1., -1., CWOBS*1.E6,   1., 1.,    1.    };
  double  CRPIX[6] =  {1.,  1., (NCHAN/2+1),     1., 1.,    1.    }; 
  double  CROTA[6] =  {0.,  0., 0.,      0., 0.,    0.,   }; 

  for(ii=0; ii<NAXIS-1; ii++){

    if(fits_make_keyn("CTYPE", ii+2, keynam, &status))
    printerror(status);
    if(fits_write_key_str(fptr, keynam, CTYPE[ii], " ", &status))
    printerror(status);

    if(fits_make_keyn("CRVAL", ii+2, keynam, &status))
    printerror(status);
    if(fits_write_key_dbl(fptr, keynam, CRVAL[ii], 10, " ", &status))
    printerror(status);

    if(fits_make_keyn("CDELT", ii+2, keynam, &status))
    printerror(status);
    if(fits_write_key_dbl(fptr, keynam, CDELT[ii],  9, " ", &status))
    printerror(status);
     
    if(fits_make_keyn("CRPIX", ii+2, keynam, &status))
    printerror(status);
    if(fits_write_key_dbl(fptr, keynam, CRPIX[ii],  9, " ", &status))
    printerror(status);

    if(fits_make_keyn("CROTA", ii+2, keynam, &status))
    printerror(status);
    if(fits_write_key_dbl(fptr, keynam, CROTA[ii],  9, " ", &status))
    printerror(status);
  }
  printf("PIXVAL etc are written....\n");

  // have to write GROUPS, GCOUNT, PCOUNT
  /* if(fits_write_key_log(fptr, "GROUPS", TRUE, " ", &status)) */
  /*   printerror(status); */
  /* if(fits_write_key_lng(fptr, "GCOUNT", GCOUNT, " ", &status)) */
  /*   printerror(status); */
  /* if(fits_write_key_lng(fptr, "PCOUNT", PCOUNT, " ", &status)) */
  /*   printerror(status); */

  //decleare ptypes
  char *PTYPE[PCOUNT] = {"UU---SIN", "VV---SIN", "WW---SIN", "BASELINE", "DATE", "DATE", "SOURCE", "FREQSEL"};
  double PSCAL[PCOUNT] = {1.e-6/FREQOBS, 1.e-6/FREQOBS, 1.e-6/FREQOBS, 1., 1., 1., 1., 1.};
  double PZERO[PCOUNT] = {0., 0., 0., 0., DATEOBS, 0., 0., 0.};
  
  printf("Now will write the PVALS....\n");

  for(ii=0; ii<PCOUNT; ii++){

    if(fits_make_keyn("PTYPE", ii+1, keynam, &status))
    printerror(status);
    if(fits_write_key_str(fptr, keynam, PTYPE[ii], "", &status))
    printerror(status);

    if(fits_make_keyn("PSCAL", ii+1, keynam, &status))
    printerror(status);
    if(fits_write_key_dbl(fptr, keynam, PSCAL[ii], 11, " ", &status))
    printerror(status);

    if(fits_make_keyn("PZERO", ii+1, keynam, &status))
    printerror(status);
    if(fits_write_key_dbl(fptr, keynam, PZERO[ii], 11, " ", &status))
    printerror(status);

 }
  printf("PSVAL etc are written....\n");
 if(fits_write_record(fptr, "         / Where baseline = 256*ant1 + ant2 + (array#-1)/10", &status))
   printerror(status);

 char HISTORY[128];

 sprintf(HISTORY, "/--------------------------------------------------------------------");
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);

 sprintf(HISTORY, "/Begin \"HISTORY\" information found in FITS tape header by UVFITS");
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);
 sprintf(HISTORY, "This data is REAL !!");
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);
 sprintf(HISTORY,"AIPS SORT ORDER = 'TB' ");
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);
 sprintf(HISTORY, "EXTEND  =                    T / FITS dataset may contain extensions");
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);
 sprintf(HISTORY, "Generated with VISFITS version        %s", versionstr);
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);
 sprintf(HISTORY, "Generated with ANTENNA FILE          %s", ANTFILE);
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);
 sprintf(HISTORY, "/End FITS tape header \"HISTORY\" information");
 if(fits_write_history(fptr, HISTORY, &status))
   printerror(status);

 if(fits_close_file(fptr, &status))
   printerror(status);
 printf("FITS header initialized...\n");
}
void Write_SU_TABLE(char *FITSFILE){

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // assumes RAAPP = SOURA and DECAPP = SOUDC
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  enum {tfields=19};
  int status ;
  fitsfile *fptr;
  char *ttype[tfields] = 
       { "ID. NO. ", "SOURCE  ", "QUAL    ", "CALCODE ", "IFLUX   ",
         "QFLUX   ", "UFLUX   ", "VFLUX   ", "FREQOFF ", "BANDWIDTH",
         "RAEPO   ", "DECEPO  ", "EPOCH   ", "RAAPP   ", "DECAPP  ",
         "LSRVEL  ", "RESTFREQ", "PMRA    ", "PMDEC   " 
       } ;
  
  char *tform[tfields] = 
    { "1J   ", "16A", "1J ", "4A ", "1E ", "1E ", "1E ", "1E ", "1D ", 
      "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D ", "1D "
    } ;
  
  char *tunit[tfields] =
    { "   ", "   ", "   ", "   ", "JY ", "JY ", "JY ", "JY ", 
      "HZ ", "HZ ", "DEGREES", "DEGREES", "YEARS  ", "DEGREES", 
      "DEGREES", "M/SEC  ", "HZ     ", "DEG/DAY", "DEG/DAY"
    } ;
 
  status=0;

  if(fits_open_file(&fptr, FITSFILE, READWRITE, &status))
    printerror(status);
  
  if(fits_create_hdu(fptr,&status))
    printerror(status);

  if(fits_write_btblhdr(fptr, NSOUR, tfields, ttype, tform, tunit, "AIPS SU", 0L, &status))
    printerror(status);

  if(fits_write_key_lng(fptr,"EXTVER", 1,"Version number of table",&status))
    printerror(status);

  if(fits_write_key_lng(fptr,"NO_IF" ,NSIDE, "The number of IFs" ,&status))
    printerror(status);
 
  if(fits_write_key_str(fptr,"VELTYP" ,"TOPOCENT", "Velocity type" ,&status)) 
    printerror(status);
  
  if(fits_write_key_str(fptr,"VELDEF" ,"RADIO", " Velocity definition" ,&status))
    printerror(status);
  
  if(fits_write_key_lng(fptr,"FREQID" ,1, "Frequency ID" ,&status))
    printerror(status);

  int ii, qual ;
  float *iflux,*qflux,*uflux,*vflux;
  double *freq_off, bandwidth;
  double *lsrvel, *rest_freq, pmra, pmdec ;
  char *source, *calc;

  iflux = (float *)calloc(NSIDE, sizeof(float));
  qflux = (float *)calloc(NSIDE, sizeof(float));
  uflux = (float *)calloc(NSIDE, sizeof(float));
  vflux = (float *)calloc(NSIDE, sizeof(float));
  
  freq_off = (double *)calloc(NSIDE, sizeof(double));
  lsrvel   = (double *)calloc(NSIDE, sizeof(double));
  rest_freq = (double *)calloc(NSIDE, sizeof(double));
  
  qual = 0; 
  bandwidth = BWOBS*1.e6;
  for(ii=0; ii<NSIDE; ii++)
    iflux[ii] = qflux[ii] = uflux[ii] = vflux[ii] = 0.;
  for(ii=0; ii<NSIDE; ii++)
    freq_off[ii] = lsrvel[ii] = rest_freq[ii] = 0.;
  pmra = pmdec = 0.;

  source = (char *)calloc(16, sizeof(char));
  calc = (char *)calloc(4, sizeof(char));

  int    SUID, row;
  char   SOUNAM[16];
  float  SOURA_HH, SOURA_Mm, SOURA_SS;
  float  SOUDC_DD, SOUDC_Mm, SOUDC_SS;
  double SOURA, SOUDC;
  double epoach;
  char   CALCOD[4], param[128];
  FILE   *fp;

  epoach = EPOACH;

  printf("Opening source file to read info....\n"); 
  fp = fopen(SOURFILE, "r");
  fgets(param,100, fp);
  fgets(param,100, fp);

  for(row=1; row<=NSOUR; row++){
    fscanf(fp, "%s%d%s%f%f%f%f%f%f", SOUNAM, &SUID, CALCOD, &SOURA_HH, &SOURA_Mm, &SOURA_SS, &SOUDC_DD, &SOUDC_Mm, &SOUDC_SS);
    
    SOURA = HMS2D(SOURA_HH, SOURA_Mm, SOURA_SS);
    SOUDC = HMS2D(SOUDC_DD, SOUDC_Mm, SOUDC_SS)/15.;
      
    if(CALCOD[0]!='C')
      sprintf(CALCOD, " ");
    sprintf(source, "%s", SOUNAM);
    strncpy(calc, CALCOD, 4);
    replace_nulls(source, 16);    
    replace_nulls(calc, 4);    
        
    if(fits_write_col(fptr, TINT,     1, row, 1,  1, &SUID,        &status))
      printerror(status);
    if(fits_write_col(fptr, TSTRING,  2, row, 1,  1,  &source,     &status))
      printerror(status);
    if(fits_write_col(fptr, TINT,     3, row, 1,  1, &qual,        &status))
      printerror(status);
    if(fits_write_col(fptr, TSTRING,  4, row, 1,  1, &calc,      &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   5, row, 1,  1, &iflux,       &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   6, row, 1,  1, &qflux,       &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   7, row, 1,  1, &uflux,       &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   8, row, 1,  1, &vflux,       &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE,  9, row, 1,  1, &freq_off,    &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 10, row, 1,  1, &bandwidth,   &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 11, row, 1,  1, &SOURA,       &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 12, row, 1,  1, &SOUDC,       &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 13, row, 1,  1, &epoach,      &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 14, row, 1,  1, &SOURA,       &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 15, row, 1,  1, &SOUDC,       &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 16, row, 1,  1, &lsrvel,      &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 17, row, 1,  1, &rest_freq,   &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 18, row, 1,  1, &pmra,        &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE, 19, row, 1,  1, &pmdec,       &status))
      printerror(status);
   }
   if(ffclos(fptr, &status))
      printerror(status);
   fclose(fp);
   printf("SU TABLE written...\n");
}
void Write_AG_TABLE(char *FITSFILE){

  enum {tfields=12};
  int status ;
  long nrows;
  double gstia0;
  fitsfile *fptr;
  char *ttype[]={ "ANNAME", "STABXYZ", "ORBPARM", "NOSTA", "MNTSTA", "STAXOF",
		  "POLTYA", "POLAA",   "POLCALA", "POLTYB","POLAB",  "POLCALB"};
  char *tform[]={"8A","3D","0D","1J","1J","1E","4A","1E","4E","4A","1E","4E"};
  char *tunit[]={" ", "METERS", " "," "," ","METERS", " ", "DEGREES", " ", " ",
		   "DEGREES", " "};
 
  status=0;
  nrows = NANTE;
  gstia0 = jd2gst(DATEOBS + REF_LON/360. + IATUTC/86400.) ;
  gstia0 = gstia0-floor(gstia0) + 1.*((int)gstia0)/360.;

  if(fits_open_file(&fptr, FITSFILE, READWRITE, &status))
    printerror(status);
  
  if(fits_create_hdu(fptr,&status))
    printerror(status);

  if(fits_write_btblhdr(fptr, nrows, tfields, ttype, tform, tunit, "AIPS AN", 0L, &status))
    printerror(status);

  if(fits_write_key_lng(fptr,"EXTVER", 1,"Version number of table",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"ARRAYX", 1657004.6290,12,"",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"ARRAYY", 5797894.3801,12,"",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"ARRAYZ", 2073303.1705,12,"",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"GSTIAO", gstia0,12,
			"GST at IAT=0 (degrees) on ref. date",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"DEGPDY",0.36098564497436661e3 ,20,
			"Earth rotation rate (deg/IAT day",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"FREQ",(FREQOBS-BWOBS/2.+CWOBS/2.)*1.E6,12,
			"Reference frequency for array",&status))
    printerror(status);
  if(fits_write_key_str(fptr,"RDATE", DATEOBS_S,
			"Reference date 'YYYY-MM-DD'",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"POLARX", 0.0,12,"Polar position X (meters) on Ref. date",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"POLARY", 0.0,12,"Polar position Y (meters) on Ref. date",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"UT1UTC", 0.0,12,"UT1-UTC (time sec.)",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"DATUTC", 0.0,12," ",&status))
    printerror(status);
  if(fits_write_key_str(fptr,"TIMSYS", "IAT",
			"Time system, 'IAT' or 'UTC' ",&status))
    printerror(status);
  if(fits_write_key_str(fptr,"ARRNAM", "GMRT","Array name",&status))
    printerror(status);
  if(fits_write_key_lng(fptr,"NUMORB", 0,
			"Number of orbital parameters",&status))
    printerror(status);
  if(fits_write_key_lng(fptr,"NOPCAL", 4,
			"Number of pol. cal constants",&status))
    printerror(status);
  if(fits_write_key_lng(fptr,"FREQID", -1 ,
  			"The ref freq of the uv data",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"IATUTC", IATUTC,12,
		       "IAT-UTC (time sec) ",&status))
    printerror(status);
  if(fits_write_key_str(fptr,"POLTYPE", "APPROX",
			"Feed polarization parameterization" ,&status))
    printerror(status);
  if(fits_write_key_lng(fptr,"P_REFANT" ,5, "Reference antenna" ,&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"P_DIFF01", 0.0,12," ???",&status))
    printerror(status);
  if(fits_write_key_dbl(fptr,"P_DIFF02", 0.0,12," ???",&status))
    printerror(status);
 
  int row;
  char  *anname, buffer1[16], buffer2[16];
  double Bxyz[3];
  long nosta,mntsta ;
  float staxof ;
  char *poltya;
  float polaa, polcala[4] ;
  char *poltyb;
  float polab, polcalb[4] ;
  FILE *fp;

  anname  = (char*)calloc(8, sizeof(char));
  poltya  = (char*)calloc(4, sizeof(char));
  poltyb  = (char*)calloc(4, sizeof(char));

  // defining common parameters
  mntsta=0;
  staxof=0.0;

  poltya[0] = 'X'; poltyb[0] = 'Y'; 
  poltya[1] = ' '; poltyb[1] = ' ' ;
  poltya[2] = ' '; poltyb[2] = ' ' ;
  poltya[3] = ' '; poltyb[3] = ' ' ;

  polaa = polab = 0.0;
  polcala[0] = polcalb[0] = 0.0;
  polcala[1] = polcalb[1] = 0.0;
  polcala[2] = polcalb[2] = 0.0;
  polcala[3] = polcalb[3] = 0.0;

  fp = fopen(ANTFILE, "r");
  for(row=1; row<=NANTE; row++){
    
    fscanf(fp, "%s%lf%lf%lf", buffer1, &Bxyz[0], &Bxyz[1], &Bxyz[2]);
    sprintf(buffer2, ":%02d", row);
    sprintf(anname, "%s%s", buffer1, buffer2);
    anname[6] = ' ';
    anname[7] = ' ';
   
    nosta = row;

    if(fits_write_col(fptr, TSTRING,  1, row, 1, 1, &anname,     &status))
      printerror(status);
    if(fits_write_col(fptr, TDOUBLE,  2, row, 1,  3, &Bxyz,    &status))
    printerror(status);
 
    //    write a null column of type 0D (0D is not ensured here
    //  if(fits_write_col_null(fptr, 3, row, 1, 1, &status))
    //  printerror(status);

 
    if(fits_write_col(fptr, TINT,     4, row, 1,  1, &nosta,        &status))
      printerror(status);
    if(fits_write_col(fptr, TINT,     5, row, 1,  1, &mntsta,        &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   6, row, 1,  1, &staxof,       &status))
      printerror(status);
    if(fits_write_col(fptr, TSTRING,  7, row, 1, 1, &poltya,     &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   8, row, 1,  1, &polaa,       &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   9, row, 1,  4, &polcala,       &status))
      printerror(status);
    if(fits_write_col(fptr, TSTRING,  10, row, 1, 1, &poltyb,     &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   11, row, 1,  1, &polab,       &status))
      printerror(status);
    if(fits_write_col(fptr, TFLOAT,   12, row, 1,  4, &polcalb,       &status))
      printerror(status);
  }
  if(ffclos(fptr, &status))
    printerror(status);
  fclose(fp);
  printf("AG TABLE written...\n");

}

void Write_FQ_TABLE(char *FITSFILE){

  int status ;
  fitsfile *fptr;
  char *ttype[] = {"FRQSEL", "IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"};
  char *tform[] = {"1J" ,    "1D" ,      "1E" ,          "1E" ,          "1J"} ;
  char *tunit[] = {"  ", "Hz", "Hz", "Hz" , "  "} ;

  status=0;
 
  if(fits_open_file(&fptr, FITSFILE, READWRITE, &status))
    printerror(status);
  
  if(fits_create_hdu(fptr,&status))
    printerror(status);

  if(fits_write_btblhdr (fptr, 1, 5, ttype, tform, tunit, "AIPS FQ", 0L, &status))
    printerror(status);

  if(fits_write_key_lng(fptr,"EXTVER", 1,"Version number of table",&status))
    printerror(status);

  if(fits_write_key_lng(fptr,"NO_IF" ,1, "The number of IFs" ,&status))
    printerror(status);

  int id, sideband;
  float cwobs, bwobs;
  double refreq;
  id = 1;
  refreq = FREQOBS*1.e6;
  refreq = 0.;
  cwobs = CWOBS*1.e6;
  bwobs = BWOBS*1.e6;
  sideband = 1;

  if(fits_write_col(fptr, TINT,    1, 1, 1, 1, &id,        &status))
    printerror(status);
  if(fits_write_col(fptr, TDOUBLE, 2, 1, 1, 1, &refreq,    &status))
    printerror(status);
  if(fits_write_col(fptr, TFLOAT,  3, 1, 1, 1, &cwobs,   &status))
    printerror(status);
  if(fits_write_col(fptr, TFLOAT,  4, 1, 1, 1, &bwobs, &status))
    printerror(status);
  if(fits_write_col(fptr, TINT,    5, 1, 1, 1, &sideband,  &status))
    printerror(status);
  
  if(ffclos(fptr, &status))
    printerror(status);
  printf("FQ TABLE written...\n");
}

void Write_scan(char *FITSFILE,  int SUID, int FREQID, int BRECORD, int ERECORD){

  long ii, indexVIS;
  int GNum;
  long pcount;
  float Time1, Time2;
  double uu, vv, ww, ha, timestamp, gst;
  float Vre, Vim, Vresky,Vimsky,gAmp, *gAmpcont, gPhi, *gPhicont, vre, vim;
  int record, status, antA, antB, gate, chan, stokes ;
  //long seedA, seedP, seedN;
  fitsfile *fptr;
  FILE *fp;
  

  pcount = PCOUNT;

  Vre = Vim = vre = vim = 0.;

  status = 0;
  if(fits_open_file(&fptr, FITSFILE, READWRITE, &status))
    printerror(status);

  //Generate constant Gain error
  gAmpcont= (float*)calloc(NTbin*NANTE, sizeof(float)); 
  gPhicont= (float*)calloc(NTbin*NANTE, sizeof(float));
 
  for(ii=0;ii<NTbin*NANTE;++ii)
    {
      gAmpcont[ii]=gsl_ran_gaussian(rA,SIG_Acont);

      gPhicont[ii]=gsl_ran_gaussian(rP,SIG_Pcont);
    }
 
  for(record=BRECORD-1; record<ERECORD; record++){
 
    /* // tjd and gst calculations has to be done...    
    timestamp = STTIM + record*INTIME/86400.;
    gst = jd2gst(DATEOBS + timestamp);
    gst = gst-floor(gst) + 1.*((int)gst)/360.;
    ha  = gst*360. + TELE_LON - souRA[SUID-1];
    ha = ha - 360.*floor(ha/360.);*/
    
    ///edit by me
    timestamp = STTIM;
    gst = jd2gst(DATEOBS + timestamp );
    gst = gst-floor(gst) + 1.*((int)gst)/360.;
    ha  = gst + TELE_LON - souRA[0];
    ha  = ha + 360.*(1.*record*INTIME)/86400.;
    ///
    
    Time1 = DATEOBS;
    //Time2 = timestamp + (1.*IATUTC)/86400.; 
    Time2 = STTIM + record*INTIME/86400. + (1.*IATUTC)/86400.; 
  
    Calc_uvw(ha, SUID);
   
    randpar[4] = 1.;
    randpar[5] = Time2;
    randpar[6] = (float)1.*SUID;
    randpar[7] = (float)1.*FREQID;
    
   /* seedA = (long) time(NULL);
    seedP = (long) time(NULL);
    seedN = (long) time(NULL);*/

    tmpbin=record*NTbin/ERECORD;
    if(tmpbin>NTbin) 
      tmpbin=NTbin;
    
    for(antA=0; antA<NANTE; antA++){
      gAmp = gAmpcont[tmpbin*NANTE+antA]+gsl_ran_gaussian(rA,SIG_A);

      gPhi = gPhicont[tmpbin*NANTE+antA]+gsl_ran_gaussian(rP,SIG_P);
      gA_re[antA] = (1.+gAmp)*cos(gPhi);
      gA_im[antA] = (1.+gAmp)*sin(gPhi);
    }

    //xshift=gsl_ran_gaussian(rPos,SIG_Pos);
    //yshift=gsl_ran_gaussian(rPos,SIG_Pos);
    //xshift=0.;
    //yshift=0.;
    for(ii=0; ii<NPSour; ii++){
      xPS[ii]=xPS0[ii]+xshift;
      yPS[ii]=yPS0[ii]+yshift;
    }
    for(antA=0; antA<NANTE; antA++)
      for(antB=antA+1; antB<NANTE; antB++){ // self are NOT written

	for(ii=0; ii<VISlen; ii++)
	  VISword[ii] = 0.;

	uu = LTT*(u[antA] - u[antB]);
	vv = LTT*(v[antA] - v[antB]);
	ww = LTT*(w[antA] - w[antB]);
	
	//random UV generator
	if(flagran==1.)
	  {
	    uu=2.*Umax*(0.5-gsl_rng_uniform(rran));
	    vv=2.*Umax*(0.5-gsl_rng_uniform(rran));
	    ww=0.;
	  }
		
	GNum = (antA+1)*256 + (antB+1);

	randpar[0] = uu;
	randpar[1] = vv;
	randpar[2] = ww;
	randpar[3] = (float)GNum;

	//write group parameters
	if(fits_write_grppar_flt(fptr, group, 1, pcount, randpar, &status ) )
	  printerror(status);

	for(chan=0; chan<NCHAN; chan++){
	  // Samirs function will come here   
	  Gen_Vis(uu, vv, ww, &Vresky, &Vimsky, chan);
	  
	  for(stokes=0; stokes<NSTOKES; stokes++){	      
	    
	      indexVIS = stokes + NSTOKES*chan;
	      
	      Vre=Vresky;
	      Vim=Vimsky;
	     
#ifdef NOISE
	      Vre += gsl_ran_gaussian(rN,SIG_N);

	      Vim += gsl_ran_gaussian(rN,SIG_N);


#endif	      
	     
#ifdef GERR
	      vre = gA_re[antA]*gA_re[antB]*Vre
		+ gA_im[antA]*gA_im[antB]*Vre
		+ gA_im[antA]*gA_re[antB]*Vim
		- gA_re[antA]*gA_im[antB]*Vim;
	      
	      vim = gA_im[antA]*gA_im[antB]*Vim
		+ gA_re[antA]*gA_re[antB]*Vim
		+ gA_re[antA]*gA_im[antB]*Vre
		- gA_im[antA]*gA_re[antB]*Vre;
	      
	      
	      VISword[RE + NCMPLX*indexVIS] = vre;
	      VISword[IM + NCMPLX*indexVIS] = vim;
#else
	      
	      VISword[RE + NCMPLX*indexVIS] = Vre;
	      VISword[IM + NCMPLX*indexVIS] = Vim;
#endif
	      VISword[WT + NCMPLX*indexVIS] = 1.;
  	    }
	}
  	

	// write visibility : one baseline
  	if(fits_write_img_flt(fptr, group, 1, VISlen, VISword, &status ))
  	  printerror(status);
	group++;
      }
  }
  if(fits_close_file(fptr, &status))
    printerror(status);
  
  printf("Total %8ld groups written till the end of scan", group-1);
 
  gsl_rng_free(rA);
  gsl_rng_free(rP);
  gsl_rng_free(rN);
  gsl_rng_free(rPos);
  //  free(VISword);
  //  free(randpar);
}
void Write_FITS_data(char *FITSFILE){

  FILE *fp;
  char param[128];
  int scan, SCNO, SUID, FREQID, BRECORD, ERECORD;

  printf("Initializing group and data array....\n");
  VISlen = NCMPLX*NSTOKES*NSIDE*NCHAN;
   
  VISword = (float*)calloc(VISlen, sizeof(float));
  randpar = (float*)calloc(PCOUNT, sizeof(float));

  gA_re = (float *)calloc(NANTE, sizeof(float));
  gA_im = (float *)calloc(NANTE, sizeof(float));

  printf("Now will write scan by scan....\n");

  fp = fopen(SCANFILE, "r");
  fgets(param,100, fp);
  fgets(param,100, fp);
  group = 1;
  for(scan=0; scan<NSCAN; scan++){
    
    fscanf(fp, "%d%d%d%d%d", &SCNO, &SUID, &FREQID, &BRECORD, &ERECORD);
    
    Write_scan(FITSFILE, SUID, FREQID, BRECORD, ERECORD);
    printf(" %2d \n", scan+1);
  }
  fclose(fp);
  GCOUNT = group-1;
  printf("FITS data written...\n");
}
void Gen_Vis(float uuc, float vvc, float wwc, float *Vre, float *Vim, int chan){
  
  int thread_id, ii, jj;
  long index;
  float thetax, thetay, ph;
  float *Vre_temp, *Vim_temp;
  float VR, VI;
  double freqcorr,nu,theta,be;

  Vre_temp = (float *)calloc(NT, sizeof(float));
  Vim_temp = (float *)calloc(NT, sizeof(float));

  freqcorr=1.+(((chan-NCHAN/2)*CWOBS)/FREQOBS);//nu/nu0
  //freqcorr=1.;//nu/nu0

  nu=FREQOBS*freqcorr*1.e6;//beam input freq. in Hz
  uuc=uuc*freqcorr;
  vvc=vvc*freqcorr;
  wwc=wwc*freqcorr;
 
  Uval=sqrt(uuc*uuc+vvc*vvc+wwc*wwc);
  VR = VI = 0.;  
  
//Diffuse part upto Uval<Umax
# ifdef DIFF
  if(Uval<Umax){
  omp_set_num_threads(NT);
# pragma omp parallel for private(thread_id,index,ii,jj,thetax,thetay,ph)
    for(index=0; index<npixels; index++){
    
      thread_id = omp_get_thread_num();
      
      ii = index%NGridI;
      jj = index/NGridI;
      
      thetax = (ii-NGridI/2)*DIx;
      thetay = (jj-NGridI/2)*DIy;
    
      ph = 2.*M_PI*(uuc*thetax + vvc*thetay+ wwc*(sqrt(1.-thetax*thetax-thetay*thetay)-1.));

      /* VR += Idata[index]*cos(ph); */
      /* VI += Idata[index]*sin(ph); */
   
      Vre_temp[thread_id] += (Idata[chan*npixels+index]*cos(ph));
      Vim_temp[thread_id] += (Idata[chan*npixels+index]*sin(ph));
    }
    
    for(ii=0; ii<NT; ii++){
      VR += Vre_temp[ii];
      VI += Vim_temp[ii];
    }
  }//check Uval<Umax loop close
#endif

#ifdef POINT

  for(ii=0; ii<NPSour; ii++){
    // xPS and yPS were in terms of DIx and DIy
    theta=sqrt((xPS[ii]*xPS[ii])+(yPS[ii]*yPS[ii]));
    be=Beam(theta,nu);
   
    ph = 2.*M_PI*(uuc*xPS[ii]+ vvc*yPS[ii] +wwc*(sqrt(1.-xPS[ii]*xPS[ii]-yPS[ii]*yPS[ii])-1.));
    VR += (IPS[ii][chan]*be*cos(ph));
    VI += (IPS[ii][chan]*be*sin(ph));
  }
#endif

  *Vre = VR;
  *Vim = VI;//phase exp(+1.*phase) 
 
  free(Vre_temp);
  free(Vim_temp);
}
int Read_Image_FITS(char *INFILE){

  fitsfile *fptr;  
  int  ii, status,  anynull,  nullval, nfound, chk;
  long fpixel, *naxes, naxis,totpixel;
  char comment[FLEN_KEYWORD];
 
  printf("Reading Image FITS file...\n");
  status = 0;
  if(fits_open_file(&fptr, INFILE, READONLY, &status))
    return status;
  if(fits_read_key_lng(fptr, "NAXIS", &naxis, comment, &status))
    return status;
  naxes = (long *)calloc(naxis, sizeof(long));
  if(fits_read_keys_lng(fptr, "NAXIS", 1, naxis, naxes, &nfound, &status))
    return status;
  
  if(fits_read_key_flt(fptr, "CDELT1", &DIx, comment, &status))
    return status;
  if(fits_read_key_flt(fptr, "CDELT2", &DIy, comment, &status))
    return status;

  if(naxes[0]!=naxes[1]){
    printf("\nWARNING: Image is not symmetric in x and y [%ld\t%ld]\n\n", naxes[0], naxes[1]); 
    return 1;
  }
  NGridI = naxes[0];
  printf("\nTotal no. of channel in  diffuse Image %ld\n",naxes[2]);
  
  printf("\nNGridI = %d\n", NGridI);

  printf("DIx    = %.2e arcmin\t", DIx);
  DIx=DR*DIx/60.;
  printf("DIx    = %.2e rad\n", DIx);
  printf("DIy    = %.2e arcmin\t", DIy);
  DIy=DR*DIy/60.;
  printf("DIx    = %.2e rad\n", DIy);

  fpixel   = 1;
  nullval  = 0;
  npixels = NGridI*NGridI;
  totpixel=npixels*naxes[2];

  Idata = (float *)calloc(totpixel, sizeof(float));
 
  if(fits_read_img(fptr, TFLOAT, fpixel, totpixel, &nullval, Idata, &anynull, &status))
    return status;
 
  if(fits_close_file(fptr, &status))
    return status;
  

  return 0;
}
int Read_Point_Sources(char *INFILE){

  FILE *fp;
  int ii,jj,ch=0;
  float tempf;

  printf("Reading Point Source file...\n");
  NPSour = 0;
  fp = fopen(INFILE, "r");
  while(!feof(fp))
    {
      ch = fgetc(fp);
      if(ch == '\n')
	NPSour++;
    }
  fclose(fp);
  
  printf("Read %d number of point sources, NCHAN=%d\n", NPSour,NCHAN);
 
  xPS = (float*)calloc(NPSour, sizeof(float));
  yPS = (float*)calloc(NPSour, sizeof(float));
  xPS0 = (float*)calloc(NPSour, sizeof(float));
  yPS0 = (float*)calloc(NPSour, sizeof(float));
  IPS = (float**)calloc(NPSour, sizeof(float*));
  IPS[0] = (float*)calloc(NPSour*NCHAN, sizeof(float));
  for(ii=1;ii<NPSour;++ii)
    IPS[ii]=IPS[0]+ii*NCHAN;
  alphapt = (float*)calloc(NPSour, sizeof(float));
  
  fp = fopen(INFILE, "r");
  for(ii=0; ii<NPSour; ii++){
    fscanf(fp, "%f%f%*f%f", &xPS0[ii], &yPS0[ii],&alphapt[ii]);//xPS[],yPS[] should be in rad.
    for(jj=0;jj<NCHAN;++jj)
      {
	fscanf(fp,"%f",&IPS[ii][jj]);
      }
  }
  fclose(fp);
 
  return 0;
}
