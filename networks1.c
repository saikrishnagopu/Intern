#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
int maxof(int a[], int num_elements);
void createsp(int w[100],int xlimit);
void interpole(double x[1400],double y[1400],double xlimit){ 
                                                               // to get the 

    int i;
        FILE *gpfp1 = fopen("interpole.gp","w");
    FILE *datafp1 = fopen("interpolo.data","w");
    if(gpfp1==NULL)
                {
            printf("Cannot generate  file\n");
                return ;
                }

   
                 for (i = 0; i<1400; i=i+1)
                {
                    
                    fprintf(datafp1, "%f %f \n",x[i], y[i]);
                   
                    }

                    fclose(datafp1);
                    fprintf(gpfp1, "set size 1,1\n");
                    fprintf(gpfp1,"set term png size 1000,700\n");
                    fprintf(gpfp1,"set border linewidth 3\n");
                    fprintf(gpfp1,"set output \"interpolediag.png\"\n");
                    fprintf(gpfp1,"plot [0:%lf]  \"interpolo.data\" \n", xlimit);
                    fclose(gpfp1);
}
void errplot(double x[1400],double y[1400],double xlimit){


    int i;
        FILE *gpfp1 = fopen("errorplt.gp","w");
    FILE *datafp1 = fopen("errorrr.data","w");
    if(gpfp1==NULL)
                {
            printf("Cannot generate  file\n");
                return ;
                }
double freq ;
   
                 for (i = 0; i<1400; i=i+10)
                {
                   freq = 1.0/x[i] ;
                    fprintf(datafp1, "%f %f \n",freq, y[i]);
                   
                    }

                    fclose(datafp1);
                    fprintf(gpfp1, "set size 1,1\n");
                    fprintf(gpfp1,"set term png size 1000,700\n");
                    fprintf(gpfp1,"set border linewidth 3\n");
                    fprintf(gpfp1,"set output \"errorgrph.png\"\n");
                    fprintf(gpfp1,"plot [0:%lf]  \"errorrr.data\" \n", xlimit);
                    fclose(gpfp1);
}
 
int
main (void)
{
  
time_t timo;
	   	double x[1400],y[1400],y1[1400],value1[2000],value2[2000],fmax;
		int i,j,k,three_pie,f[100],w[100],levels;
		int wmax;
        printf("Enter the levels : ");
			three_pie = 10 ; 
			srand((unsigned) time(&timo));
		 
			for(i=0;i<5;i++)
				 w[i] = rand() % 15 ;
			scanf("%d",&levels);	 
			createsp(w,4);
			value1[0] = 0 ;
	    	value2[0] = 0 ;
			for(i=1;i<levels/2;i++)
			 value1[i] = value1[i-1] + (10.0/levels) ;
			for(i=1;i<levels/2;i++)
			 value2[i] = value2[i-1] - (10.0/levels) ;

			wmax = maxof(w,5);
			fmax = wmax/6.5 ;
			 x[0]=0.0; 
printf("Max ferquency is : %f",fmax);
double tmin;
tmin = 1.0/(2.1*fmax);
	  
			for(i=0;i<1400;i++){ 
		    k=0; 
			y[i] = sin(w[0]*x[i])+sin(w[1]*x[i])+sin(w[2]*x[i])+sin(w[3]*x[i])+sin(w[4]*x[i]);
		    y1[i] = y[i];
			if(y[i]>0){
			while(y[i]>value1[k])
			k++;		 
			y[i] = value1[k-1]; 
			  }
		else if(y[i]==0){
		y[i] = 0.0;
	}
		   else{        
				while(y[i] < value2[k])
				 k++;
			 y[i] = value2[k-1];
			 }

			x[i+1] = x[i] + tmin ;
}
    
double err,err1,ms1error;

err = 10.0 / levels ;
int l;
srand((unsigned) time(&timo));
	i=0;
  y[0]=0.0;
	while(i<1400){
   
     err1 = err *( (double)rand()/RAND_MAX );
     err1 = err1 * (((double)rand()/RAND_MAX)-((double)rand()/RAND_MAX));
	 y[i+50] = y[i] + err1 ;     
	 i= i+50;
   
	}
     
int N;
N=1399;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  const gsl_interp_type *t = gsl_interp_cspline_periodic; 
  gsl_spline *spline = gsl_spline_alloc (t, N);

  double xi, yi,oreo,mserror;
  double xero[1400],yero[1400],msero[1400];
  printf ("#m=1,S=0\n");
  gsl_spline_init (spline, x, y, N);
xi=0.0;
  for (i = 0; i <= 1400; i++)
    {
      
      xi = xi+tmin/150.0;
      
      yi = gsl_spline_eval (spline, xi, acc);
      xero[i]=xi;
      yero[i] = yi;
      oreo = sin(w[0]*x[i])+sin(w[1]*x[i])+sin(w[2]*x[i])+sin(w[3]*x[i])+sin(w[4]*x[i]);
      mserror += pow(oreo-yero[i],2);
      msero[i] = mserror/(xi) ;
    }

mserror = mserror/ (x[1399]-x[0]);
                   mserror = sqrt(mserror);
  printf("mean square error is : %f\n",mserror);
yero[0]=0;
 interpole(xero,yero,4.0);
errplot(xero,msero,100.0);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return 0;
}

		int maxof(int a[100], int num_elements)
		{
		   int i, max=0;
		   for (i=0; i<num_elements; i++)
		   {
			 if (a[i]>max)
			 {
				max=a[i];
			 }
		   }
		   return(max);
		}

		void createsp(int w[100],int xlimit){
		int i=0;
		FILE *gpfp = fopen("suppos.gp","w");
		fprintf(gpfp, "set size 1,1\n");
		fprintf(gpfp,"set term png size 1000,700\n");
		fprintf(gpfp,"set border linewidth 3\n");
		fprintf(gpfp,"set output \"suposition.png\"\n");
		fprintf(gpfp,"plot [0:%d] ", xlimit);

		for(i=0;i<5;i++)
		{
		fprintf(gpfp,"sin(%d*x)",w[i]);
		if(i<4) fprintf(gpfp,"+");}
		fclose(gpfp);

		}

