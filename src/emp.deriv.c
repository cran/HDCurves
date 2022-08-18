#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/*****************************************************************************************
 * The following are the inputs of the function that are read from R
 *
 * Input:
 * nobs = vector whose entries indicate number of observations per subject
 * y = vector containing nobs values
 * t = vector containing nobs time (time is incremented by one) values for each player
 * ztype = integer indicating what sequence of quotient differences to use
    0 - Mine (employing both D and u)
    1 - Cote (employing only u)
    2 - De Brabantar.
 * u = double for the u parmater in the empirical derivative formulation
 * D = integer for the D paramter in the empirical derivative formulation
 *
 * Ouput:
 * z = vector of size nobs where empirical derivative values are stored
 *****************************************************************************************/


void empderiv(int *nobs, double *y, double *t, int *ztype, int *D,  double *u, double *z){

  emp_derivative(y, t, *D, *u, *nobs, *ztype, z);

}
