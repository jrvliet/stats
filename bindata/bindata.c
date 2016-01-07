
// Command line arguements:
// $1 = input file name
// $2 = column data are in
// $3 = input data type (0=linear, 1=log10)
// $4 = bin size (equal log10)
// $5 = lower bin limit (log10)
// $6 = upper bin limit (log10)


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[]) {


    int maxdata = 4000000;
    int i, j, col, ndata, nbins, type, nhdr;
    double zbrent;
    double n, CL, xllim, xullim, lgxcen, lgxlo, lgxhi;
    double x[maxdata], bin[maxdata];
    double binup[maxdata], bindn[maxdata];
    double binsize, z1, z2, lo, hi, xx;
    char in_file[200], out_file[200];
    char string[200], newline[200];
    double iondense, size, columnDense;
    double kpctocm = 3.086e21; 
    double column[maxdata];
    double dum;

    // Assume a 1-sigma single sided confidence level
    CL = 0.8413;

    // Assume a tolerance of 1e-5 for numerical root solving for
    // the Poisson error; the maximum value of any bin is XULIM;
    // if we exceed this we will need to catch it and terminate
    tol = 1.0e-7;

    // Read in command line arguements
    strcpy(in_file, argv[1])
    col = (int)argv[2]
    type = (int)argv[3]
//    binsize = (double)argv[4]
//    z1 = (double)argv[5]
//    z2 = (double)argv[6]

    type = 0;
    z1 = 10.0;
    z2 = 24.0;
    binsize = 0.5;

    // Cell size column = 9
    // Ion density column =14

    // Build the outfile
    strcpy(out_file, in_file);
    strcat(out_file, ".logfreqbin");

    // Read in the data file, store in x[i]
    FILE *fp = fopen(in_file, "r");
    fgets(newline, sizeof(newline), fp);
    i = 0;
    ndata = 0;
    while(fgets(newline, sizeof(newline),fp)){

        sscanf(newline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &dum,  &dum,  &dum,  &dum,  &dum,  &dum,  &dum,  &dum,  &size, 
                &dum,  &dum,  &dum,  &dum,  &dum,  &iondense);
        columnDense = size * kpctocm * iondense;
        column[i] = columnDense;
        ndata += 1;
    }
    fclose(fp);

    // Bin the data
    // Datum is in a a bin if it is >= lower bound and < the 
    // upper bound, but don't miss xmax; store normilzation 
    // possibilities on the fly (MAXBIN, NORM); compute the
    // uncertainty in the bin value (use Poisson stats)
    
    // set counters
    nbins = (int)( (z2-z1)/binsize);

    // Loop over bins
    for( i=0; i<nbins; i++){

        // Initialize
        n = 0;
        bin[0] = 0;

        lo = pow(10.0, z1) * pow(10.0, (i-1)*binsize);
        hi = pow(10.0, z1) * pow(10.0, i*binsize);

        // Loop over data and count how many are in the bins
        for( j=0; j<ndata; j++){
                        
            data = column[j];
            if (data>=lo && data<hi){
                bin[i]+=1;
            }
            if (i==nbins-1 && data==hi){
                bin[i] += 1;
            }
            bindn[i] = 0;
            binup[i] = 0;
        }

        // Compute Poisson errors
        // If bin is large, then use Gaussian approx
        if (bin[i]>200.0){
            binup[i] = bin[i] + sqrt(bin[i]);
            bindn[i] = bin[i] - sqrt(bin[i]);
        }
        else{
            if (bin[i]>0 ){
                n = bin[i]
                xllim = 0.0;
                xulim = 1.0e3*n;
                bindn[i] = zbrent(CLdn, xllim, xulim, tol);
                binup[i] = zbrent(CLup, xllim, xulim, tol);
            }
        }
    }


    // Normalize and output
    FILE *fpout = fopen(out_file, 'w');
    
    fprintf(fpout, "Input data file %s\n", in_file);
    fprintf(fpout, "Bin center\tlog(F)\t|bin/2|\tdlogFdn\tdlogFup\n");
    // Loop over the bins
    for (i=0; i<nbins; i++){  
    
        lo = pow(10.0, z1) * pow(10.0, (i-1)*binsize);
        hi = pow(10.0, z1) * pow(10.0, i*binsize);

        lgxlo = log10(lo);
        lgxhi = log10(hi);
        lgxcen = 0.5*(lgxhi+lgxlo);
        
        // Frequency distribution normiization (per x veriable per Ndata)
        bin[i] = bin[i] / (hi-lo) / (float)(ndata);
        bindn[i] = bindn[i] / (hi-lo) / (float)(ndata);
        binup[i] = binup[i] / (hi-lo) / (float)(ndata);

        // Convert errors to log10
        if (bin[i]>0){
            bindn[i] = log10(bindn[i]/bin[i]);
            binup[i] = log10(binup[i]/bin[i]);
        }
        else{
            bindn[i] = 0;
            binup[i] = 0;
        }

        fprintf(fpout, "%lf \t %lf \t %lf \t %lf \t %lf \n",
                    lgxcen, bin[i], 0.5*(lgxhi-lgxlo), bindn[i], binup[i]);

    }
    
    fclose(fpout);

    return 0;
}




























