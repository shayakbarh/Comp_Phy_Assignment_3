#include<stdio.h>
#include<math.h>
#include<gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include <complex.h>




double sinc(double x) {
    if (x == 0.0) {
        return 1.0;
    } else {
        return sin(x) / x;
    }
}

// Function for plotting two curves using Gnuplot
void plotTwoCurves(double *x1, double *y1, double *x2, double *y2, int size) {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe) {
        // commands to Gnuplot to plot two curves
        fprintf(gnuplotPipe, "plot '-' with lines title 'Numerical FT', '-' with lines title 'Analytic FT'\n");
        
        //  data for the first curve
        for (int i = 0; i < size; ++i) {
            fprintf(gnuplotPipe, "%lf %lf\n", x1[i], y1[i]);
        }
        fprintf(gnuplotPipe, "e\n");
        
        // data for the second curve
        for (int i = 0; i < size; ++i) {
            fprintf(gnuplotPipe, "%lf %lf\n", x2[i], y2[i]);
        }
        fprintf(gnuplotPipe, "e\n");

        fflush(gnuplotPipe);
        printf("Press enter to exit...\n");
        getchar(); // Wait for user to press enter

        // Close the pipe
        pclose(gnuplotPipe);
    } else {
        printf("Error opening Gnuplot pipe.\n");
    }
}

// Function for analytical fourier transfom
double* analytical_FT(double *x, int size) {
    double* y = malloc(size * sizeof(double));
    if (y == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < size; ++i) {
        if (x[i] >= -1.0 && x[i] <= 1.0) {
	    y[i] = sqrt(M_PI/2.0);
        } else {
            y[i] = 0.0;
        }
    }
    return y;
}



// Function to roll an array of doubles by half of its length and return the rolled array
double* roll_array(double arr[], int size) {
    // Allocate memory for the rolled array
    double *rolled_array = (double*) malloc(sizeof(double) * size);

    // Find the mid point
    int mid = size / 2;

    // Roll the array
    for (int i = 0; i < size; i++) {
        rolled_array[i] = arr[(i + mid) % size];
    }

    return rolled_array;
}

// Function to perform element-wise complex multiplication and exponentiation
double complex* calculate_FT(double *karr, double *xarr, double *dft_arr, int size) {
    // Allocate memory for the result array
    double complex *result = (double complex*) malloc(sizeof(double complex) * size);

    // Perform element-wise multiplication, multiplication by i, and complex exponentiation
    for (int i = 0; i < size; i++) {
        double complex temp = -I*karr[i]*xarr[0];
        result[i] = (xarr[1]-xarr[0]) * sqrt(size/(2*(M_PI))) *cpow(M_E, temp)* ((dft_arr[2*i])+(I*dft_arr[2*i+1]) ) ;
    }

    return result;
}

// Function to calculate 2*pi*x/(n*dx) for each element in the array
double* calculate_k(double *x, int n) {
    // Allocate memory for the result array
    double *result = (double*) malloc(sizeof(double) * n);
    
    // Calculate expression for each element
    for (int i = 0; i < n; i++) {
        
      result[i] = 2 * M_PI *(i -n/2) / (n * (x[1]-x[0]));
      
          
    }

    // Return the array containing the calculated values
    return roll_array(result,n);
}



int main()
{
int i;
double x;
 int n=256;
 double x_values[n];
 double y_values[n];
double data[2*n];
double x_min=-50;
double x_max=50;
 double FT[n];
double delta=(x_max-x_min)/(n-1);
for(i=0;i<n;i++)
{
  
x=x_min+i*delta;
 x_values[i]=x;
 y_values[i]=sinc(x);
 REAL(data,i)=sinc(x);
IMAG(data,i)=0.0;
}

gsl_fft_complex_radix2_forward(data,1,n);

 for(i=0;i<n;i++)
{
  data[2*i]= data[2*i]/sqrt(n);
  data[2*i+1]=data[2*i+1]/sqrt(n);
}

    // Calculate k values
    double *k_values = calculate_k(x_values,n);
    // Calculate FT values
    complex *FT_values = calculate_FT(k_values,x_values,data,n);

    // Calculate analytic values
    double *f_analytic = analytical_FT(k_values,n);

    
    
    for (int i = 0; i < n; i++) {
      printf(" %.5lf\n", creal(FT_values[i]) );
    }
    

    for (int i = 0; i < n; i++) {
      FT[i]=creal(FT_values[i]) ;
	
    }
    
    // Open a file for writing

    FILE *file = fopen("gsl_FT_sinc.txt", "w");

    if (file == NULL) {

        fprintf(stderr, "Error opening file.\n");

        return 1;

    }



    // Print the results to the file

    for (int i = 0; i < n; i++) {

        fprintf(file, "%f\n", creal(FT_values[i]));

    }



    // Close the file


    
    //plot the two curves
    plotTwoCurves(k_values,FT,k_values,f_analytic,n);

    return 0;
 
return 0;
}
