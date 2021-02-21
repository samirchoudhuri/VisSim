//Generate diffuse image at different freq. in FITS fromat
// uses input angular power spectrum of brightness temperature fluctuation 

#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<fitsio.h>
#include<unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

fftw_complex *in;
double *out;

int N,NB;
double length;

double CRVAL[3],CDELT[3],CRPIX[3];
long naxes[3],naxesim[3];
void SWRITE_HDR(char *);

double Beam(double theta,double freq);
void power_spec(double *,double *,int *);

double A_150,alpha=2.8,betav,nu0;
double k_B=1.38e3 /* in Jy*/,c=3.0e8;


//return power spectrum
double  P_I(double nu,double U)//freq. in Hz
{
  double pu;
  pu=pow(2.*k_B*nu*nu/(c*c),2.)*A_150*pow(nu0/nu,2.0*alpha)*pow(1000./(2.*M_PI*U),betav);
  //pu=A_150*pow(1000./(2.*pi*U),betav);
  
  /////////////////////////////////////////////////////////////

  //pu=pow(2.*k_B*nu*nu/(c*c),2.); // unit angular power spectrum 

  /////////////////////////////////////////////////////////////

  return(pu);
}

int main(int argc, char *argv[])
{
  fftw_plan p,p1;
  int *num;
  double *binu,*ps_nobeam,*ps_beam,*img;
  double p11=0.,p22=0.;
  char OPT[1],psnbeam[150],pswbeam[150];
  int i,j,k,index,index1,xdim,ydim,ia,Nchan;
  double L,pixel,nu,deltanu,chan0,fac,u,amp,spindex,delspindex;
  double thetax,thetay,theta;
  double theta0,A;
  double be;

  //for random number generator
  unsigned long int seed;
  gsl_rng *r;
  double sigma=1.;
  //done

  FILE *fp,*fp1;
  
  // reads name of parameter file and output FITS file at runtime
  if(argc!=3)
   {
      printf("Usage: %s <input> <ouput FITS file>\n", argv[0]);
      return 1;
    }

  if(access(argv[2], F_OK)==0)
    {
      printf("File %s  exists\n",argv[2]);
      exit(0);
    }
 
  //reading input parameter
  // seed, No. of grid points, resolution in arcmin, central freq., channel reolution., chan0, number of channel, 
  // No. of bin in power spec, 
  //theta0 for power spectrum normalization, sp. index for freq. scaling A_150 betav
  //Name of power spectrum data file without and with PB

  fp=fopen(argv[1],"r");
  fscanf(fp,"%ld%d%lf%lf%lf%lf%d%d%lf%lf%lf%lf%lf",&seed,&N,&L,&nu0,&deltanu,&chan0,&Nchan,&NB,&theta0,&spindex,&delspindex,&A_150,&betav);
  fscanf(fp,"%s%s",psnbeam,pswbeam);
  fclose(fp);
  
  printf("N=%d nu0=%e deltanu=%e chan0=%lf Nchan=%d\nNB=%d theta0=%e\tspindex=%lf delspindex=%lf seed=%ld\n",N,nu0,deltanu,chan0,Nchan,NB,theta0,spindex,delspindex,seed); 
  printf("%s %s %s betav=%lf\n",pswbeam,psnbeam,argv[2],betav);
  ydim=(N/2+1);
  xdim=N;

  pixel=L;// pixel resolution in arcmin
  L=M_PI*L/(180.*60.);// in rad
  theta0=M_PI*theta0/(180.*60.);// theta0 in rad
  A=M_PI*theta0*theta0/2.;//Normalization constant for power spectrum
  fac=L/(sqrt(2.)*N);length=L*N;// factor for sp. intensity fluctuations
  //fac=1./(sqrt(2.)*N*L);
  printf("pixel=%earcmin (%erad)\nA=%e fac=%e length=%e\n",pixel,L,A,fac,length);
 

  //memory allocation for out and in  array
  out=(double*)calloc ((N*(N+2)),sizeof(double));
  in=(fftw_complex*)&out[0];
  img=(double*)calloc(N*N,sizeof(double));
 
  p= fftw_plan_dft_c2r_2d (N, N, in,out, FFTW_ESTIMATE);
  p1= fftw_plan_dft_r2c_2d (N, N, out,in, FFTW_ESTIMATE);

  num=(int*)calloc(NB,sizeof(int));
  binu=(double*)calloc(NB,sizeof(double));
  ps_nobeam=(double*)calloc(NB,sizeof(double));
  ps_beam=(double*)calloc(NB,sizeof(double));

  //for random number generators
  r= gsl_rng_alloc(gsl_rng_cmrg);
  gsl_rng_set (r, seed);
  //done

  //Filling Fourier Components  
  //along axis (j-0 and j=N/2)
  for(j=0;j<ydim;j=j+N/2)
    for(i=1;i<N/2;++i)
      {
	// along + x 
	u=sqrt(1.*(i*i+j*j))/length;
	amp=fac*sqrt(P_I(nu0,u));
	index=i*ydim+j;
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	
	// along -x 
	index1=(N-i)*ydim+j;
	in[index1][0]=in[index][0];
	in[index1][1]=-in[index][1];
	
      }
  // upper half plane excluding x axis
  
  for(i=0;i<xdim;++i)
    for(j=1;j<N/2;++j)
      {
	ia= (i>N/2) ? (N-i) : i ;
	u=sqrt(1.*(ia*ia+j*j))/length;
	amp=fac*sqrt(P_I(nu0,u));
	index=i*ydim+j;
	
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
      }
      
  //4 points remain 
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      {
	if(i+j==0) 
	  {
	    in[0][0]=0.0;
	    in[0][1]=0.0;
	  }
	else
	  {
	    u=(N/2.)*sqrt(1.*(i*i+j*j))/length;
	    amp=fac*sqrt(P_I(nu0,u));
	    index=i*(N/2)*ydim+j*(N/2);
	    
	    in[index][0]=pow(-1.,(i*N/2+j*N/2))*amp*gsl_ran_gaussian(r,sigma);
	    in[index][1]=0.0;
	  }
      }
  
  // finished filling Fourier components

  //power spectrum estimate without PB
  power_spec(binu,ps_nobeam,num);
  
  //FFT of complex arry
  fftw_execute(p);
  
  //Writing in FITS format
  //Header information
  printf("Creating new file and writing header\n\n");
  naxesim[0]=(long) N;naxesim[1]=(long) N;naxesim[2]=(long) Nchan;
  CRVAL[0]=0.;CRVAL[1]=0.;CRVAL[2]=nu0;
  CDELT[0]=1.*pixel;CDELT[1]=1.*pixel;CDELT[2]=deltanu;
  CRPIX[0]=(1.*N)/2.+1;CRPIX[1]=(1.*N)/2.+1;CRPIX[2]=chan0;
  

  //Write FITS header
  SWRITE_HDR(argv[2]);
  
  fitsfile *fptrim;
  long fimpixel[3],limpixel[3];
  int status=0;


  //Write multichannel data in FITS format with freq. scaling (PB and Sp. intensity) 
  fits_open_file(&fptrim,argv[2],READWRITE,&status);
  for(k=0;k<Nchan;k++)
    {
      nu=nu0+(k+1.+0.5-chan0)*deltanu;
      fimpixel[0]=1;fimpixel[1]=1;fimpixel[2]=(long) (1+k);
      limpixel[0]=(long)N;limpixel[1]=(long)N;limpixel[2]=(long)(1+k);
      
      printf("chan=%d nu=%e image=[%ld %ld %ld] [%ld %ld %ld]\r",k,nu,fimpixel[0],fimpixel[1],fimpixel[2],limpixel[0],limpixel[1],limpixel[2]);
      for(i=0;i<N;++i)
	for(j=0;j<N;++j)
	  {
	    index=i*(N+2)+j;
	    index1=j*N+i;
	    thetax= (i-N/2)*L;
	    thetay= (j-N/2)*L;
	    theta=sqrt((thetax*thetax)+(thetay*thetay));
	    be=Beam(theta,nu);
	    //be=Beam(theta,nu0);
	    //be=1.;
	    img[index1]=be*out[index]*pow((nu0/nu),(spindex+(delspindex*log10(nu/nu0))));//to generate Power spec using beam multiplied out[] array
	  }
      fits_write_subset(fptrim,TDOUBLE,fimpixel,limpixel,img,&status);
    }

  //To calculate power spectrum after multiplied with beam
  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      {
	index=i*(N+2)+j;
	thetax= (i-N/2)*L;
	thetay= (j-N/2)*L;
	theta=sqrt((thetax*thetax)+(thetay*thetay));
	be=Beam(theta,nu0);
	out[index]=be*out[index];//to generate Power spec using beam multiplied out[] array
      }
 
  //FFT of real array 
  fftw_execute(p1);

  //calculate Power Spec after multiplying with PB
  power_spec(binu,ps_beam,num);
  

  //Write power spectrum in data file
  fp=fopen(psnbeam,"w");
  fp1=fopen(pswbeam,"w");
  for(i=0;i<NB;++i)
    if(num[i]>0)
      {
  	fprintf(fp,"%e %d ",binu[i],num[i]);
  	fprintf(fp1,"%e %d ",binu[i],num[i]);
	p11=ps_nobeam[i]*pow(N,4)/(length*length);
	fprintf(fp,"%e %e\n",P_I(nu0,binu[i]),p11);
	p22=ps_beam[i]/A;//for 610 cont=7.22e-5
	fprintf(fp1,"%e %e\n",P_I(nu0,binu[i]),p22);
      }
  fclose(fp);
  fclose(fp1);
  
  gsl_rng_free(r);
  free(binu);
  free(num);
  free(ps_beam);
  free(ps_nobeam);
  fits_close_file(fptrim,&status);
  fftw_destroy_plan(p);
  fftw_destroy_plan(p1);
  fftw_free(out);
}


//Function to calculate the power spectrum
void power_spec(double *usum,double *pksum,int *no)
{
  double Umax,Umin,binsiz,u;
  int i,j,index,ia,tmp;

  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
    }

  Umax=sqrt(2.)*(N/2.)/length;
  Umin=1./(length);
  binsiz=(log10(1.*N/sqrt(2.))/(1.*NB));
  
  for(i=0;i<N;++i)
    for(j=0;j<N/2+1;++j)
      {
  	index=i*(N/2+1)+j;
  	ia=(i>N/2) ? (N-i) : i ;
  	u=sqrt(1.*(ia*ia+j*j))/length;
  	if((u>Umin)&&(u<Umax))
  	  {
  	    tmp=(int) floor(log10(u/Umin)/binsiz);
	    tmp=(tmp<NB) ? tmp: tmp-1;
  	    no[tmp]++;
	    pksum[tmp]+=((in[index][0]*in[index][0])+(in[index][1]*in[index][1]));
  	    usum[tmp]+=u;
  	  }
      }
  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]=usum[i]/no[i];
	pksum[i]=pksum[i]/no[i];
      }
}
