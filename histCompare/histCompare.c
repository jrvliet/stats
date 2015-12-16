
// 
// Filename: histCompare.c
// Author: Jacob Vander Vliet
// Date: 2/20/15
//
// Reads in histograms of the phase of absorbing cells
// Performs the chi^2 test from Numerical Recipies in C (Sec. 14.3)
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <vicNl.h>

#define ITMAX 1000       // Maximum allowed number of iterations.
#define EPS 3.0e-7      // Relative accuracy.
#define FPMIN 1.0e-30   // Number near the smallest representable floating-point number.
#define DOLOG 0
#define NORM 0         // Set to zero always

/*
Given the arrays bins1[1..nbins] and bins2[1..nbins] , containing two sets of binned data, and given the number of constraints knstrn (normally 1 or 0), this routine returns the number of degrees of freedom df, the chi-square chsq, and the significance prob. A small value of prob indicates a significant difference between the distributions bins1 and bins2. Note that bins1 and bins2 are both float arrays, although they will normally contain integer values.
*/

void chstwo(float bins1[], float bins2[], int nbins, int knstrn, float *df, float *chsq, float *prob) {

  float gammq(float a, float x);
  int j;
  float temp;
  *df=nbins-knstrn;
  *chsq=0.0;
  
  float sum1 = 0.;
  float sum2 = 0.;
 
  for (j=1; j<=nbins; j++){
    sum1 += bins1[j];
    sum2 += bins2[j];
  }

  for (j=1;j<=nbins;j++){
    if (bins1[j] == 0.0 && bins2[j] == 0.0)
      --(*df);
    else {
      temp = sqrt(sum1/sum2)*bins2[j] - sqrt(sum2/sum1)*bins1[j];
      *chsq += temp*temp/(bins1[j]+bins2[j]);
    }
  }
  *prob=gammq(0.5*(*df),0.5*(*chsq));
}



/*
Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
*/
float gammq(float a, float x){
  void gcf(float *gammcf, float a, float x, float *gln);
  void gser(float *gamser, float a, float x, float *gln);
  void nrerror(char error_text[]);
  float gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}


/*
Returns the incomplete gamma function P (a, x) evaluated by its series representation as gamser. Also returns ln Γ(a) as gln.
*/
void gser(float *gamser, float a, float x, float *gln) {
  float gammln(float xx);
  void nrerror(char error_text[]);
  int n;
  float sum,del,ap;

  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } 
  else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    nrerror("a too large, ITMAZ too small in routine gser");
    return;
  }
}

/*
Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction represen-
tation as gammcf. Also returns ln Γ(a) as gln.
*/

void gcf(float *gammcf, float a, float x, float *gln) {
  float gammln(float xx);
  void nrerror(char error_text[]);
  int i;
  float an,b,c,d,del,h;
  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

float gammln(float xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}

static char vcid[] = "$Id: nrerror.c,v 3.2 1999/08/23 23:59:06 vicadmin Exp $";

/* Numerical Recipes standard error handler */
void nrerror(char error_text[])
{
	void _exit();
	fprintf(stderr,"Model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	_exit(1);
}



int main() {

  int numBins = 2500;
  int numIons = 4;
  int numGals = 3;
  int i, j, k, l;
  int knstrn = 0;

  float ionMean, ionMeans[numIons];
  float runMeans[3];
  float df, chsq, prob;
  float number, sum;
  float hists[numGals][numBins];


  char filename[30];
  //  char *fileloc = "/home/matrix3/jrvander/sebass_gals/dwarfs/abscells/";  
  char *fileloc = "/home/jacob/matrix/sebass_gals/dwarfs/abscells/";
  char file0[100];
  char line[50];
  
  FILE *fp;

  // List of ions
  char *ion_list[numIons];
  ion_list[0] = "HI";
  ion_list[1] = "MgII";
  ion_list[2] = "CIV";
  ion_list[3] = "OVI";

  // List of galaxy IDs
  char *galID_list[numGals];
  galID_list[0] = "D9o2";
  galID_list[1] = "D9q";
  galID_list[2] = "D9m4a";
  
  // List of runs
  char *run_list[numGals];
  run_list[0] = "SN";
  run_list[1] = "ALL_1";
  run_list[2] = "ALL_8";
  
  // Loop thru ions
  for (i=0; i<numIons; i++){

    // Loop thru gals
    for (j=0; j<numGals; j++){

      // Construct the filename
      strcpy(filename, galID_list[j]);
      strcat(filename, ".");
      strcat(filename, ion_list[i]);
      strcat(filename, ".abscell.hist");
      strcpy(file0, fileloc);
      strcat(file0, filename);
   
      // Read in the histograms
      fp = fopen(file0, "r");
      k = 0;
      while (fscanf(fp, "%f", &number) != EOF){
	if (DOLOG == 1){
	  if (number>0.){
	    hists[j][k] = log10(number);
	  }
	  else {
	    hists[j][k] = 1e-10;
	  }
	}
	else {
	  hists[j][k] = number;
	}
	k++;  
      }
      fclose(fp);

    }

    // If NORM is high, then normalize the histograms
    if (NORM==1){
      for (k=0; k<numGals; k++){
	sum = 0.;
	for (l=0; l<numBins; l++){
	  sum += hists[k][l];
	}
	
	// Normalize 
	for (l=0; l<numBins; l++){
	  hists[k][l] /= sum;
	}
      }
    }
	
	



    ionMean = 0.0;

    printf("Ion: %s\n", ion_list[i]);

    chstwo(hists[0], hists[1], numBins, knstrn, &df, &chsq, &prob);
    printf("Comparing %s to %s:\n", run_list[0], run_list[1]);
    printf("Chi^2:              %f\n", chsq);
    printf("Degrees of Freedom: %f\n", df);
    printf("Probability:        %f\n", prob);
    printf("Reduced Chi^2:      %f\n\n", chsq/df);
    ionMean += chsq/df;
    
    chstwo(hists[0], hists[2], numBins, knstrn, &df, &chsq, &prob);
    printf("Comparing %s to %s:\n", run_list[0], run_list[2]);
    printf("Chi^2:              %f\n", chsq);
    printf("Degrees of Freedom: %f\n", df);
    printf("Probability:        %f\n", prob);
    printf("Reduced Chi^2:      %f\n\n", chsq/df);
    ionMean += chsq/df;

    chstwo(hists[1], hists[2], numBins, knstrn, &df, &chsq, &prob);
    printf("Comparing %s to %s:\n", run_list[1], run_list[2]);
    printf("Chi^2:              %f\n", chsq);
    printf("Degrees of Freedom: %f\n", df);
    printf("Probability:        %f\n", prob);
    printf("Reduced Chi^2:      %f\n\n", chsq/df);
    ionMean += chsq/df;    

    /*
    chstwo(hists[1], hists[1], numBins, knstrn, &df, &chsq, &prob);
    printf("Comparing %s to %s:\n", run_list[1], run_list[1]);
    printf("Chi^2:              %f\n", chsq);
    printf("Degrees of Freedom: %f\n", df);
    printf("Probability:        %f\n", prob);
    printf("Reduced Chi^2:      %f\n\n", chsq/df);
    */

    printf("\n");
    ionMean /= 3.;
    ionMeans[i] = ionMean;

  }


  printf("\n\n");
  for (l=0; l<numIons; l++){
    printf("Chi^2 mean for %s: \t %f\n", ion_list[l], ionMeans[l]);
  }


  return 0;
}
