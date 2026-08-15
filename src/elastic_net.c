/**
 * Filename: elastic_net.c
 * 
 * Description: Implements elastic net regression fit similar to glmnet 
 * (Driver for functions in elastic_net_functions.c)
 * 
 * Version 2 - Added warm start optimization
 * Version 3 - Added edits for active set optimization
 * Version 4 - Added covariance updates
 * Version 5 - Cleaned up comments, added references to thesis equations
 * 
 * Author: Marc Lane
 * Version: 5.0
 * Date: 2026-08-12
 */

#include <stdlib.h>
#include <string.h>
#include "elastic_net_functions.h"
#include "elastic_net.h"

/* Elastic-net regularisation for a fixed alpha, using covariance updates.
 *
 * Inputs:
 *   Z_in          - n x p column-major covariate matrix
 *   z_in          - length-n target station vector
 *   n             - number of days
 *   p             - number of covariate stations
 *   lambdas       - list lambda values options from CV
 *   n_lambda      - number of lambda values
 *   alpha         - elastic-net mixing hyperparameter
 *   thresh        - convergence threshold, compared against the squared change in beta
 *   maxit         - cap on the total number of sweeps
 *   beta_out      - n_lambda x p, gets filled with the fitted coefficients
 *   intercept_out - length-n_lambda, gets filled with the intercept for each lambda
 *   C             - p x p covariance matrix, allocated by ENCE_impute.c or MICE_impute.c
 *   has_row       - length-p flags for which rows of C have been computed, allocated by ENCE_impute.c or MICE_impute.c
 *
 *   Returns 0 on success, otherwise returns:
 *     1 - for a memory allocation error
 *     2 - for invalid inputs
 *     3 - if maxit was reached by coord_descent_cov()
 */
int elastic_net_cov(const double *Z_in, const double *z_in, int n, int p, const double *lambdas, int n_lambda, double alpha,
                          double thresh, int maxit, double *beta_out, double *intercept_out, double *C, char *has_row)
{

    /* Return 2 for invalid inputs */
    if (n <= 0 || p <= 0 || n_lambda <= 0 || thresh <= 0.0 || maxit <= 0 || alpha < 0.0 || alpha > 1.0 || !Z_in ||
        !z_in || !lambdas || !beta_out || !intercept_out || !C || !has_row) {
        return 2;
    }

    /* Allocating memory for buffers and copies */
    double *Z    = malloc((size_t)n * p * sizeof(double));
    double *z    = malloc((size_t)n * sizeof(double));
    double *Zm   = malloc((size_t)p * sizeof(double));
    double *Zs   = malloc((size_t)p * sizeof(double));
    double *beta = calloc((size_t)p, sizeof(double));  // calloc to set beta_k's to zero at the beginning
    int *is_active  = malloc((size_t)p * sizeof(int));
    int *active_idx = malloc((size_t)p * sizeof(int));
    double *g = malloc((size_t)p * sizeof(double));

    /* Allocation error check */
    if (!Z || !z || !Zm || !Zs || !beta || !is_active || !active_idx || !g) {
        free(Z); free(z); free(Zm); free(Zs); free(beta); free(is_active); free(active_idx); free(g);
        return 1;
    }

    /* Making copies of data */
    memcpy(Z, Z_in, (size_t)n * p * sizeof(double));
    memcpy(z, z_in, (size_t)n * sizeof(double));

    /* Perform standardisation transformation once for all lambdas*/
    double zm, zs;
    if (standardise(Z, z, n, p, Zm, Zs, &zm, &zs) != 0) {
        free(Z); free(z); free(Zm); free(Zs); free(beta); free(is_active); free(active_idx); free(g);
        return 2; // return 2 if standardise function fails (due to invalid inputs)
    }

    /* g = c when all beta = 0 at the beginning (Equation 2.3.10) */
    for (int k = 0; k < p; ++k) {

        const double *col_k = Z + (size_t)k * n;  // pointer to the start of the k-th column of Z

        /* Initialise g_k = c_k for all k */
        double c_k = 0.0;
        #pragma omp simd reduction(+:c_k) // vectorisation reduction pragma to speedup this dot product
        for (int i = 0; i < n; ++i){
            c_k += col_k[i] * z[i];   // Equation 2.3.11
        }
        g[k] = c_k; // Equation 2.3.10 reduces to this when beta=0 at the start of the algorithm
    }


    /* Loop each lambda in the CV lambda options */
    for (int l = 0; l < n_lambda; ++l) {

        /* beta, g and C carry over from the previous lambda */
        int cd_status = coord_descent_cov(Z, beta, n, p, lambdas[l] / zs, alpha, thresh, maxit, is_active, active_idx, g, C, has_row);

        /* coord_descent_cov() returns 1 if maxit was reached */
        if (cd_status == 1) {
            free(Z); free(z); free(Zm); free(Zs); free(beta); free(is_active); free(active_idx); free(g);
            return 3;
        }
        
        /* coord_descent_cov() returns 2 if invalid inputs were inputted */
        if (cd_status == 2) {
            free(Z); free(z); free(Zm); free(Zs); free(beta); free(is_active); free(active_idx); free(g);
            return 2;
        }
        

        /* Back-transform a copy of beta so beta stays standardised for the next lambda's warm start */
        double *bcol = beta_out + (size_t)l * p;  // pointer to the start of this lambda's column of beta_out
        memcpy(bcol, beta, (size_t)p * sizeof(double)); // copy each beta from each coord_descent_cov for each lambda into bcol

        /* Back-transform beta to the original scale */
        intercept_out[l] = back_transform(bcol, Zm, Zs, zm, zs, p);
    }

    free(Z); free(z); free(Zm); free(Zs); free(beta); free(is_active); free(active_idx); free(g);
    
    return 0;

}
