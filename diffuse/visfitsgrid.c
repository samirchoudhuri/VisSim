// This program reads a multichannel FITS image (diffuse) and UVFITS template then, add gridded visibility to the nearest uv track
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>
# include <fftw3.h>
# include "read_fits_func.h"
# include "fitsprog.h"

extern double CRVAL[],CDELT[],CRPIX[];
extern long naxes[],naxesim[];

extern float del_chan,nu_chan0,chan0,Umaxc;
extern long nstokes,nchan,ncmplx,gcount;

int main(int argc, char *argv[])
{ 
  char IMAGEFITS[128],UVFITS[128];
  int chan,n_ave,ii,jj,ii1,jj1;
  int nbasln,N;
  long index,index1;
  double dU,fac;
  
 
  if(argc!=4)
    {
      printf("Usage: %s <input FITS image file> <input UVFITS file> <factor=0(overwrite, 1=add)>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[1],"%s",IMAGEFITS);
  sscanf(argv[2],"%s",UVFITS);
  fac=atof(argv[3]);

  if(access(IMAGEFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",IMAGEFITS);
    exit(0);
  }

  if(access(UVFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",UVFITS);
    exit(0);
  }

  printf("Image Header\n");
  SREAD_HDR(IMAGEFITS);
  N=naxesim[0];n_ave=naxesim[2];
  dU=(double) ((180.*60.)/(M_PI*N*CDELT[0]));

  printf("pixel=%lf arcmin naxes image={%ld,%ld,%ld}\ndU=%e\n",CDELT[0],naxesim[0],naxesim[1],naxesim[2],dU);
  printf("N=%d n_ave=%d\n",N,n_ave);

  //read image and FFT and fill in a array
  fitsfile *fptr;
  long fpixel[naxis],lpixel[naxis],inc[naxis];
  int status=0,anynull;
  double nulval=0;

  fftw_complex *in,*reim;
  double *out,*img;
  fftw_plan p1;

  reim=(fftw_complex*)calloc ((N*(N/2+1)*n_ave),sizeof(fftw_complex));
  out=(double*)calloc ((N*(N+2)),sizeof(double));
  in=(fftw_complex*)&out[0];
  img=(double*)calloc(N*N,sizeof(double));
  p1= fftw_plan_dft_r2c_2d (N, N, out,in, FFTW_ESTIMATE);

  printf("\n");

  fits_open_file(&fptr,IMAGEFITS,READONLY,&status);
  for(chan=0;chan<n_ave;++chan)
    {
      fpixel[0]=1;fpixel[1]=1;fpixel[2]=(long)(chan+1);
      lpixel[0]=N;lpixel[1]=N;lpixel[2]=(long)(chan+1);
      inc[0]=1;inc[1]=1;inc[2]=1;;
      printf("image read  [%ld,%ld,%ld] [%ld,%ld,%ld]\r",fpixel[0],fpixel[1],fpixel[2],lpixel[0],lpixel[1],lpixel[2]);
      
      fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,&nulval,img,&anynull,&status);
      
      for(jj=0;jj<N;++jj)
  	for(ii=0;ii<N;++ii)
  	  {
  	    index=jj*N+ii;
  	    index1=ii*(N+2)+jj;
  	    out[index1]=img[index];
  	  }
     
      fftw_execute(p1);
      for(ii=0;ii<N;++ii)
  	for(jj=0;jj<=N/2;++jj)
  	  {
  	    index1=chan*(N/2+1)*N+ii*(N/2+1)+jj;
  	    index=ii*(N/2+1)+jj;
  	    reim[index1][0]=pow(-1.,ii+jj)*in[index][0];//Image centre at (N/2,N/2)
  	    reim[index1][1]=pow(-1.,ii+jj)*(-1.*in[index][1]);//r2c use -1, to make it +1 exponent multiply with -1
  	  }
    }
  fftw_destroy_plan(p1);
  fits_close_file(fptr,&status);

  printf("\n");

  //Only to read UVFITS header(switch on to write rms values,ff,Umaxc in read_fits_func.c)
  nbasln=readfits(UVFITS,0.,0.,1.,1.);
 
  //put grid value in UV track
  float randpar[3],*data;
  long group;
  long nel,el1=1;
  double uuc,vvc,signv;
  int stokes,anynul=0;
  float nulval1=0;
  double *corrfact;

  nel=ncmplx*nstokes*nchan;
  printf("nel=%ld gcount=%ld\n",nel,gcount);
  
  printf("n_ave=%d (image) nchan=%ld (vis)\n",n_ave,nchan);
  if(n_ave!=nchan)
    {
      printf("n_ave should be equal to nchan\n");
      return 1;
    }
  
  corrfact = (double*) calloc(n_ave, sizeof(double));
  for(ii=0; ii<n_ave; ii++)
    {    
      corrfact[ii] = 1.+ (del_chan/nu_chan0)*(ii+1.+0.5-chan0);//lamda[chan0]/lambda[ii]
    }
 
  data=(float*)calloc(nel,sizeof(float));

  fits_open_file(&fptr,UVFITS,READWRITE,&status);
  for(group=1;group<=gcount;group++)
    {
      fits_read_grppar_flt(fptr,group,el1,3,randpar,&status);
      //rand parameter read done    
      signv= (randpar[1]<0.) ? -1. : 1. ;
            
      fits_read_img_flt(fptr,group,el1,nel,nulval1,data,&anynul,&status);
      
      for(chan=0;chan<n_ave;chan++)
	for(stokes=0;stokes<nstokes;++stokes)
	  {
	    index=chan*nstokes+stokes;
	    data[3*index+0]*=fac;
	    data[3*index+1]*=fac;
	  }
      
      for(chan=0;chan<n_ave;chan++)
	{
	  uuc=signv*randpar[0]*corrfact[chan];
	  vvc=signv*randpar[1]*corrfact[chan];
	  
	  //uuc=signv*randpar[0]*corrfact[chan]*nu_chan0;
	  //vvc=signv*randpar[0]*corrfact[chan]*nu_chan0;

	  if(abs(uuc)<(N*dU/2.) && abs(vvc)<(N*dU/2.))
	    {
	      ii1 = (int)roundf(uuc/dU);
	      ii=(ii1<0) ? N+ii1 : ii1 ;
	      jj = (int)roundf(vvc/dU);
	      
  	      index1=chan*(N/2+1)*N+ii*(N/2+1)+jj;
	      
  	      for(stokes=0;stokes<nstokes;++stokes)
  		{
  		  index=chan*nstokes+stokes;
  		  data[3*index+0]+=reim[index1][0];
  		  data[3*index+1]+=(signv*reim[index1][1]);
		}
  	    }
  	}
      fits_write_img_flt(fptr, group, 1, nel, data, &status );
    }
  if(fits_close_file(fptr, &status))
    printerror(status);
}
