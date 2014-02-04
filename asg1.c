#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#define pi 3.1456
#define divisions 10000

double Sigval (double Am[], double Fr[], double x){
	double signal=0;
	int i;

	for (i=0; i<5; ++i){
	signal+=(Am[i]*sin(2* pi*Fr[i]*x));}

	return signal;
}

void quantise(double s[],double xi[],double Fr[],double Am[],double ts,int c){
int i=0;
for(i=0;i<c;i++){ s[i]= (double)(floor((Sigval(Am,Fr,i*ts))*100)/100);
                  xi[i] = i*ts;}


}


void distort(double s[],double mina ,int c){

    int i=0;
    for(i=0;i<c;i++){
            s[i]+= pow(-1,rand()%2)*((double)rand()/(double)(RAND_MAX))*mina;// printf("%g\n",s[i]);

                    }
}

void interpolation(double y[],double x[],int cnt,double e[],double Am[],double Fr[],double ts){

{
  double xi,yi;
FILE *f;
f =fopen("dataout.dat","w");
if(f==NULL){printf("file not found\n");}


    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline  = gsl_spline_alloc (gsl_interp_cspline, cnt);
gsl_spline_init (spline, x, y, cnt);

    int i=0;
    for (xi = 0.0; xi <= x[cnt-1]; xi += 0.01) {
//printf("Entered here: %g--%g\n",xi,x[cnt-1]);
        yi = gsl_spline_eval (spline, xi, acc);
       e[i++]=  pow((Sigval(Am,Fr,xi)-yi),2);
        fprintf (f,"%g		%g\n", xi, yi);
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
fclose(f);
  }

}


double msr(double e[],int pcnt){
  int i=0;double er=0.0;
for(i=0;i<pcnt;i++){er+= e[i];}
 return er/pcnt;
}



int main()
{
    double fmax,Amax,fs,x; double Am[5],Fr[5];int i=0;

    printf("Enter frequency max:");
    scanf("%lf",&fmax);
    //Amax=fmax;
  printf("Enter Amplitude max:");
    scanf("%lf",&Amax);

    srand((unsigned int)time(NULL));
    for(i=0;i<5;i++){ Fr[i] = ((double)rand()/(double)(RAND_MAX))*fmax ; //printf("Fr[%d]= %g",i,Fr[i]);
                       Am[i] = ((double)rand()/(double)(RAND_MAX))*Amax ; //printf("Am[%d]= %g \n",i,Am[i]);
                    }


double ts, MaxFr=0, MinFr=10000,rms = 1000;
	double TPeriod;
	double MaxAm=0, MinAm = 10000;




	for (i= 0; i< 5 ; i++){
	if (Fr[i] > MaxFr)	MaxFr =  Fr[i];
	if (Fr[i] < MinFr)	MinFr =  Fr[i];
	}
   fs= 2*MinFr;
  // int cnt= (int)MaxFr/fs;
int cnt =30;

    double s[cnt],xi[cnt];
	TPeriod = 1/MinFr;
       ts = 1/fs;


 quantise(s,xi,Fr,Am,ts,cnt);
     FILE *f;int pcnt=0;
f= fopen("data.dat","w");
if(f==NULL){printf("file not found\n");}
for(x=0.0;x<=xi[cnt-1];x=x+0.01){
//printf("came here\n");
pcnt++;
fprintf(f,"%g		%g\n",x,Sigval(Am,Fr,x));
}
fclose(f);



	double p;double e[pcnt];
	for(i=0; i< divisions ; i++){
	x = i*TPeriod/divisions;
	p= Sigval ( Am, Fr, x);
	p = abs(p);
	if (p>MaxAm) MaxAm = p;
	if (p<MinAm) MinAm = p;
	}

	//L = ceil(MaxAm);
       // printf("levels are : %d\n",L);
	double mina=100000;
	for (i=0; i< 5;++i){
	if (mina > Am[i]) mina= Am[i];
	}
	//printf("Minimum Amplitude %g\n",mina);
     mina=0.01;

f= fopen("errdata.dat","w");
if(f==NULL){printf("file not found\n");}

     for(fs=2*MinFr;  fs< 4*fmax ;fs+=0.01){
            //printf("camr ehre\n");
     ts = 1/fs;
    quantise(s,xi,Fr,Am,ts,cnt);
     distort(s,mina,cnt);
  interpolation(s,xi,cnt,e,Am,Fr,ts);
   rms= msr(e,pcnt); 
   fprintf(f,"%f        %f\n",fs,rms);
              }
              printf("errror we got is: %f\n",rms);
fclose(f);

printf("the value of fmax is obtained by seeing the graph\n ");
  return 0;
}

