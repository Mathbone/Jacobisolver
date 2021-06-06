 /*
 ========================================================================

 This program is an iterative solver of Poisson equation
 u_xx + u_yy = f in [x0,x1]x[y0,y1];

 ========================================================================
 input:
 (defined in the programm)
         Nx,Ny          - the number of grid points in x- and y-directions;
         x0,x1,y0,y1    - the boundaries of the computational domain in x-
                          y-directions;
         
 output: a,b            - are the Jacobi iteration approximations to the 
                          solution of the Poisson equation;
 
 ========================================================================
 
Brandon Wilson, 12/07/2017, ISU, Pocatello, ID
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/time.h>



void rhs_bc_anltsol(int Nx, int Ny, double hx, double hy,
                    double *a, double *b, double *UE);

double normL2_R(int N, double *x);
double cpuSecond(void) ;


int main(int argc, char **argv){
    
    int l, Nx, Ny, K;
    double x0,y0,x1,y1,hx,hy,norm;
    double *UE, *a, *b;
    double *res;
	double *next, *last, *tmp;
    clock_t time0, time1;
    double start_t_n, end_t_n, compute_t_n;

/*
     Computational Domain
*/
    
    x0 = 0;
    y0 = 0;
    
    x1 = 1;
    y1 = 1;
    
    
    /*
     Nx, Ny - the number of grid points in x- and y- directions,
     K = Nx*Ny - the size of the unknown vector, hx, hy - the grid
     sizes in x- and y- directions.
     */
    
    Nx=600;
    Ny=Nx;
    
    printf("\n \tThe size of the matrix is %d by %d \n", Nx, Ny );
    
    K=(Nx+2)*(Ny+2);
    hx=(x1-x0)/(Nx+1);
    hy=(y1-y0)/(Ny+1);
    
    
    /*
     Allocate and define analytic solution UE, and initialize 
     arrays a and b with the initial guess.
     */
    
    UE  = malloc(K * sizeof(double));
    a = malloc(K * sizeof(double));
    b = malloc(K * sizeof(double));
	next = a;
	last = b;
    
    rhs_bc_anltsol(Nx,Ny,hx,hy,a,b,UE);

    res  = malloc(K * sizeof(double));
    
    double coef = 1/(2/(hy*hy)+2/(hx*hx));
	double hxx = 1/hx/hx;
	double hyy = 1/hy/hy;
	double maxerr = 0;
        
    time0 = clock();
    start_t_n = cpuSecond();
	int count =0;
    
    do{
		for(int j = 1; j<Ny+1; j++){
			int my = j*(Nx+2);
			int mym = (j-1)*(Nx+2);
			int myp = (j+1)*(Nx+2);
			for(int i = 1; i < Nx+1; i++){
				next[i+my] = coef*(-1+hxx*(last[i-1+my]+
						last[i+1+my])+hyy*(last[i+mym]+last[i+myp]));
			}
		}
		maxerr =0;
		for(l = 0; l < K; l++){
			if(fabs(a[l]-b[l])>maxerr){
				maxerr = fabs(a[l]-b[l]);
			}	
		}

		tmp = next;
		next = last;
		last = tmp;
		count++;
	}while(maxerr>=.0001 && count<100000);
    
    end_t_n = cpuSecond();
    time1 = clock();
    compute_t_n = end_t_n - start_t_n;
    
    printf("\n\tSolver time = %f sec \n",(float)(time1-time0)/CLOCKS_PER_SEC);
    printf("   \tSolver wall-time (sys/time) = %f sec \n\n", compute_t_n);
    

    
    /*
     Calculate the l2-norm of the residual
    */

    
     printf("\titerations =  %d \n\n",count);

    /*
     Calculate the error between the approximate and analytic solution a-UE
    */
    
    
    for (l=0; l < K; l++){
        res[l] = last[l] - UE[l];
    }

     norm = normL2_R(K, res);
     printf("\t||error||_2 =  %10.7e \n\n",norm);
    
    
    free(UE);
    UE = NULL;
    free(a);
    a = NULL;
	free(b);
	b = NULL;
    

    free(res);
    res = NULL;    	     	 	 
    
    
    return 0;
} /* end function main */


double cpuSecond(){
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ( (double)tp.tv_sec + (double) tp.tv_usec*1.e-6 );
}


