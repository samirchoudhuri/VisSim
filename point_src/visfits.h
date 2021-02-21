// Definitions
# define VERSION "0.01"  // format A.AA 
# define IATUTC (34.)

# define REF_LON 82.5
# define C 299792458.
# define DR (M_PI/180.)
# define EPOACH (2000.)

//static double sqrarg;
//#define SQR(a) ((sqrarg=(a)) == 0.D ? 0.D : sqrarg*sqrarg)
static double cubarg;
#define CUB(a) ((cubarg=(a)) == 0.D ? 0.D : cubarg*cubarg*cubarg)

// CONSTANTS GEneral
enum{RE=0, IM=1, WT=2};


//CONSTANTS for FITS
enum{NANTE=30, NSTOKES=2, NCMPLX=3, PCOUNT=8, NAXIS=7, NSIDE=1};

// funcs.c functions
void get_version();
void Legals();
void printerror(int);
void Calc_uvw(double, int);
void Read_Inputs(char *);
void Init_PARM();
void Init_HDR(char *);
void Write_SU_TABLE(char *);
void Write_FQ_TABLE(char *);
void Write_AG_TABLE(char *);
void Write_scan(char *, int, int, int, int);
void Write_FITS_data(char *);
int Read_Image_FITS(char *);
void Gen_Vis(float, float, float, float *, float *,int );
int Read_Point_Sources(char *);

// utils.c functions
void CHKFILE(char filename[]);
void replace_nulls(char *, int);
double jd2gst(double );
char *jd2iau_date(double);
double HMS2D(float, float, float);
void COPYLEFT();


# define MAX_VAL (10000000.)

# define COPYLEFTFILE "COPYLEFT.INFO"


/* You can comment a particulr thing here to get read of that effect */
/* For example if you do " // #define POINT " then the point source  */
/* part of the code will not be compiled */

//# define POINT
//# define DIFF
# define NOISE
//# define GERR
