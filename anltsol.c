#include <stdlib.h>
#include <math.h>

/* #define M_PI acos(-1.0) */

void rhs_bc_anltsol(int Nx, int Ny, double hx, double hy,
                   double *a, double *b, double *UE)
{

/* Calculate the r.h.s for the test problem and the corresponding
   analytic solution UE.
========================================================================
 input:  Nx,Ny - the number of grid points in x- and y-directions.
         hx,hy - the grid steps in x- and y-directions,
 
 output: rhs      - the right hand side of the discretized system,
         UE       - the analytic solution at the grid points,
========================================================================

 Dr. Yury Gryazin, 03/01/2017, ISU, Pocatello, ID
*/


int i,j,my;
double x,y;

/*
  The analytic solution and bou
*/


    for( j = 0; j < Ny+2; j++) {
        y=(j)*hy;
        my = j*(Nx+2);

        for( i = 0; i < Nx+2; i++){
            x = (i)*hx;
            
            UE[i + my]  = .25*( x*(x-1) + y*(y-1) );
			if(i == 0 || i== Nx+1  || j == 0 || j ==Ny+1){
            	a[i + my]  = .25*( x*(x-1) + y*(y-1) );
            	b[i + my]  = .25*( x*(x-1) + y*(y-1) );
			}
			else{
				a[i+my] = 0;
				b[i+my] = 0;
			}
        }
    }

    
    
    
return;
}
