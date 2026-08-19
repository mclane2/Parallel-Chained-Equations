/**
 * Filename: elastic_net_functions.c
 * 
 * Description: This file contains the functions called in elastic_net.c
 * 
 * Equations numbers reference my thesis: "Chained Equations for Incomplete Climate Networks"
 * 
 * Function descriptions:
 * update_beta()       - Updates beta_k, called only by coord_descent_cov()
 * update_g()          - Updates g, called only by coord_descent_cov()
 * coord_descent_cov() - Main cyclical coordinate descent function
 * standardise()       - Performs standardisation transformation on matrix Z and target vector z
 * back_transform()    - Transforms standardised scale coefficients to original scale coefficients
 * 
 * Version 2 - Added active set cycling
 * Version 3 - Added SIMD vectorisation pragmas
 * Version 4 - Added covariance updates
 * Version 5 - Cleaned up comments, added references to thesis equations
 * 
 * Warm-starting carries beta, g, C and has_row across the lambda values in CV rather than being reset every fit. 
 * 
 * 
 * Author: Marc Lane
 * Version: 5.0 
 * Date: 2026-08-19
 */
#include <stdlib.h>
#include <math.h>

#include "elastic_net_functions.h"


/* Beta_k Update
 *
 * Given g_k and the current beta_k, computes the new beta_k and returns the change
 *
 * Inputs:
 *   beta_k - pointer to the current value of beta_k
 *   g_k    - current value of g_k
 *   l1     - l1 (= lambda * alpha ) is precomputed in coord_descent_cov() to avoid it recomputing it everytime function is called
 *   l2     - l2 (= lambda * (1 - alpha)) is also precomputed for the same reason
 * 
 */
static inline double update_beta(double *beta_k, double g_k, double l1, double l2)
{

    double beta_old = *beta_k;
    /* Equation 2.3.8 */
    double u = g_k + beta_old;

    /* Soft-threshold update (Equation 2.3.3 and 2.3.4) */
    if (u > 0.0 && l1 < fabs(u)) {
        *beta_k = (u - l1) / (1.0 + l2);
    } else if (u < 0.0 && l1 < fabs(u)) {
        *beta_k = (u + l1) / (1.0 + l2);
    } else {
        *beta_k = 0.0;
    }

    /* Return change in Beta */
    return *beta_k - beta_old;
}


/* Updates g after beta_k changes by delta using: g = g - delta*C_k
 *
 * If station k has never been active before, compute row k of C
 *
 * Inputs: 
 *   Z     - n x p standardised matrix of covariate stations (Z_{-j})
 *   n     - number of days
 *   p     - number of covariate stations
 *   k     - index of the current station
 *   delta - change in beta
 *   g     - lenght-p, vector of g_k values
 *   C     - p x p covariance matrix
 *   has_row  - length-p vector indictating which rows of C have already been computed
 * 
 */
static inline void update_g(const double *Z, int n, int p, int k, double delta, double *g, double *C, char *has_row)
{
    double *C_k = C + (size_t)k * p; // pointer to the start of the k-th row of C 

    /* if this row of C has not been computed before */
    if (!has_row[k]) {

        const double *col_k = Z + (size_t)k * n;  // pointer to the start of the k-th column of Z

        /* Compute row k of C as the dot product of column k with every other column, including itself */
        /* Loop through the other stations using index l to compute the dot products */
        for (int l = 0; l < p; ++l) {

            const double *col_l = Z + (size_t)l * n;  // pointer to the start of the l-th column of Z

            /* Computing entry C_kl of C*/
            double C_kl = 0.0;
            #pragma omp simd reduction(+:C_kl) // vectorisation reduction pragma to speedup this expensive dot product
            for (int i = 0; i < n; ++i){ 
                C_kl += col_k[i] * col_l[i]; // Computing the dot product C_kl (Equation 2.3.12)
            }
            C_k[l] = C_kl; // Fill in the l-th entry of the k-th row of C
        }

        /* This row of C has been computed, flagging it for reuse */
        has_row[k] = 1;
    }

    /* Using the k-th row of C (either just computed or getting reused), update g */
    #pragma omp simd // vectorisation pragma for speeding up this computation
    for (int l = 0; l < p; ++l){
        g[l] -= delta * C_k[l];  // Equation 2.3.14
    }
}


/* Main cyclical coordinate descent function for elastic-net regression, using covariance updates
 *
 * Minimises the elastic net objective function for standardised data
 *
 * Inputs:
 *   Z          - n x p column-major covariate matrix
 *   beta       - length-p model coefficient vector
 *   n          - number of days
 *   p          - number of covariate stations
 *   lambda     - regularisation strength hyperparameter
 *   alpha      - elastic-net mixing hyperparameter
 *   threshold  - convergence threshold determined by the maximum value of delta^2
 *   maxit      - maximum value of the inner+outer loop iterations (glmnet default = 100,000)
 *   is_active  - length-p, active-set membership flags for each station (1=active station, 0=non-active station)
 *   active_idx - length-p, indices of each of the members of the active set
 *   g          - length-p, current g_k for every station. 
 *   C          - p x p, covariance matrix, with each entry equal to C_k,l = <Z_k, Z_l>
 *   has_row    - length-p vector indictating which rows of C have already been computed
 *
 * Returns 0 if converged, 1 if maxit was hit, 2 if invalid inputs were given
 */
int coord_descent_cov(const double *Z, double *beta, int n, int p, double lambda, double alpha, double threshold, int maxit, 
                      int *is_active, int *active_idx, double *g, double *C, char *has_row)
{

    /* Return 2 early for invalid inputs */
    if (n <= 0 || p <= 0 || maxit <= 0 || lambda < 0.0 || alpha < 0.0 || alpha > 1.0) {
        return 2;
    }

    /* Precompute l1 and l2 */
    const double l1 = lambda * alpha;
    const double l2 = lambda * (1.0 - alpha);

    /* Active set starts empty */
    int n_active = 0;
    for (int k = 0; k < p; ++k){
        is_active[k] = 0;
    }
    
    int sweeps = 0;     // Total sweeps so far

    /* Outer loop, sweeps until betas stop changing */
    double dlx = threshold;
    while (dlx >= threshold && sweeps < maxit) {

        dlx = 0.0;      // Convergence tracker for full sweep

        /* For each covariate station k (the other stations) */
        for (int k = 0; k < p; ++k) {

            /* Update beta using Equation 2.3.3, returns delta = new_beta_k - old_beta_k */
            double delta = update_beta(&beta[k], g[k], l1, l2);

            /* If this beta_k changed*/
            if (delta != 0.0) {

                /* Update g to reflect the change in beta */
                update_g(Z, n, p, k, delta, g, C, has_row);

                /* Update the largest delta^2 (dlx) if delta^2>dlx */
                if (delta * delta > dlx){
                    dlx = delta * delta;
                }
            }

            /* If this Beta is nonzero and not in the active set, then add it to the active set */
            if (beta[k] != 0.0 && !is_active[k]) {
                is_active[k] = 1;
                active_idx[n_active] = k;
                n_active++;     
            }
        }

        sweeps++; // Count this full sweep

        /* new delta^2 max value tracker for the active set loop */
        double dlx_a = threshold;

        /* Inner loop, sweep only the active set until it converges */
        while (dlx_a >= threshold && sweeps < maxit) {

            dlx_a = 0.0;

            /* For each active covariate station */
            for (int a = 0; a < n_active; ++a) {

                const int k = active_idx[a];    // Index of the a-th active station

                /* Update beta with 2.3.3, returns delta = new_beta_k - old_beta_k */
                double delta = update_beta(&beta[k], g[k], l1, l2);

                /* If this beta_k changed */
                if (delta != 0.0) {

                    /* Update g to reflect the change in beta */
                    update_g(Z, n, p, k, delta, g, C, has_row);

                    /* Update delta^2 for this sweep */
                    if (delta * delta > dlx_a){
                        dlx_a = delta * delta;
                    }
                }
            }

            ++sweeps;   // Count this active sweep
        }
    }

    /* Return 0 if converged, return 1 if maxit was reached */
    return (sweeps >= maxit) ? 1 : 0;
}


/* Standardises Z and z
 *
 * Set each column of Z and the target vector z to mean zero, then rescales each so that its sum of squares is 1 (Equation 2.3.1). 
 * The original means and standard deviations are returned so that back_transform() can reverse this afterwards.
 *
 * This function assumes no column of Z and no target z is constant, since a constant vector has zero std dev and scaling it would divide by zero.
 * This should be checked before this function is called.
 *
 * Inputs:
 *   Z  - n x p matrix of covariate stations (Z_{-j})
 *   z  - length-n target station vector (Z_j)
 *   n  - number of days
 *   p  - number of covariate stations
 *   Zm - length-p, gets filled with the original mean of each column of Z
 *   Zs - length-p, gets filled with the original std dev of each column of Z
 *   zm - gets filled with the original mean of z
 *   zs - gets filled with the original std dev of z
 *
 * After call:
 *   Every column of Z has sum 0 and sum of squares 1
 *   z has sum 0 and sum of squares 1
 * 
 * Returns 1 for invalid inputs, returns 0 for success
 * 
 */
int standardise(double *Z, double *z, int n, int p, double *Zm, double *Zs, double *zm, double *zs)
{
    /* Return 1 early for invalid inputs */
    if (n <= 0 || p <= 0) {
        return 1;
    }

    /* Looping through each column */
    for (int j = 0; j < p; ++j) {

        double *col = Z + (size_t)j * n; // pointer to the j-th column of Z

        /* Computing the mean of the j-th column */
        double sum = 0.0;
        for (int i = 0; i < n; ++i){
            sum += col[i];
        }
        Zm[j] = sum / n;

        /* Subtract mean, accumulate sum of squared deviations */
        double ssq = 0.0;
        for (int i = 0; i < n; ++i) {
            col[i] -= Zm[j];
            ssq += col[i] * col[i];
        }

        /* Computing std dev */
        Zs[j] = sqrt(ssq / n);

        /* Scale so column has sum of squares 1 (Equation 2.3.1) */
        const double scale = sqrt(ssq);
        for (int i = 0; i < n; ++i){
            col[i] /= scale;
        }
    }

    /* Computing mean of the target vector z */
    double sum = 0.0;
    for (int i = 0; i < n; ++i){
        sum += z[i];
    }
    *zm = sum / n;

    /* Subtract mean, accumulate sum of squared deviations */
    double ssq = 0.0;
    for (int i = 0; i < n; ++i) {
        z[i] -= *zm;
        ssq += z[i] * z[i];
    }

    /* Computing std dev */
    *zs = sqrt(ssq / n);

    /* Scale so z vector has sum of squares 1 (Equation 2.3.1) */
    const double z_scale = sqrt(ssq);
    for (int i = 0; i < n; ++i){
        z[i] /= z_scale;
    }

    return 0;
}


/*
 * Converts standardised coefficients back to original-scale coefficients and returns the intercept.
 *
 * The solver fits on centred data, so it needs no intercept. 
 * Undoing the centring reintroduces one, equal to zm minus the sum of Zm[j]*beta[j] over every station.
 *
 * Inputs:
 *   beta - fitted coeffficients on the standardised scale
 *   Zm   - column means of original Z   (from standardise)
 *   Zs   - column std devs of original Z (from standardise)
 *   zm   - mean of original z
 *   zs   - std dev of original z
 *   p    - number of covariate stations
 *
 * After Call:
 *   beta transformed to original scale coefficients
 *   Returns regression intercept
 */
double back_transform(double *beta, const double *Zm, const double *Zs, double zm, double zs, int p)
{
    double a = zm;

    for (int j = 0; j < p; ++j) {
        beta[j] *= zs / Zs[j];      // back-transform to original scale
        a -= Zm[j] * beta[j];       // accumulating the intercept, using the rescaled beta[j]
    }
    return a;
}
