double simp(double (*func)(double),int n,double a,double b)
{       
        int i;
        double y=0.0;
        double sump,sumq,sumr,h;

        h = (b-a)/(1.0*n); 
        sump = func(a)+func(b);
	sumq=0.0;
	sumr=0.0;
	for(i=1;i<n;i+=2)
           {          
           sumq=sumq+func(a+i*h);
            }
	for(i=2;i<n-1;i+=2)
          {
	    sumr=sumr+func(a+i*h);
	  }
	y=  h*(sump + 4*sumq+ 2*sumr)/3.0 ;
      
      return(y);

}

