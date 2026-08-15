/**
 * Filename: MICE_DURR_impute.c
 * 
 * Description: Implements the elastic_net_DURR and MICE_get_params R functions in C to allow for OpenMP parallelism
 * 
 * Function descriptions:
 * predict_mse()       - Computes the mean-squared-error (MSE) for determining the best performing hyperparameter in CV
 * column_impute()     - Imputes all missing values of a target station by fitting an elastic net model
 * mice_impute_cycle() - A single chained equations cycle over each target station (with synchronous updates)
 * 
 * 
 * Version 1 - Edited from ENCE_impute.c for MICE DURR
 * Version 2 - Added covariance updates
 * Version 3 - Cleaned up comments, added references to thesis equations
 * 
 * Author: Marc Lane
 * Version: 3.0
 * Date: 2026-08-15
 */
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "elastic_net.h"
#include "MICE_DURR_impute.h"



/* Computes the mean squared error of the elastic-net prediction
 *
 * Inputs:
 *   X         - n x p covariate matrix
 *   beta      - coefficients vector
 *   intercept - model intercept
 *   y         - length-n observed targets
 *
 * Returns MSE
 */
static double predict_mse(const double *X, int n, int p, const double *beta, double intercept, const double *y)
{
    double sse = 0.0;
    /* Loop through each row of X */
    for (int i = 0; i < n; ++i) {

        /* For each entry of beta, compute the prediction using Equation 2.2.1 */
        double yhat = intercept;
        for (int j = 0; j < p; ++j) {
            yhat += X[(size_t)j * n + i] * beta[j];
        }
        /* The error is the true value minus the prediction */
        double e = y[i] - yhat; 
        sse += e * e;
    }
    return sse / n;
}


/* Imputes all missing values of the target station by fitting an elastic net model
 *
 *   df_old           - n x p data matrix
 *   target_out       - length-n, gets filled with output column for this target station
 *   n                - number of days
 *   p                - number of stations
 *   col              - index of the target station
 *   missing_col      - length-n, nonzero values were originally missing for this target station
 *   fold_col         - length-n, CV fold id for observed days 
 *   lambdas/alphas   - CV hyperparameter options list
 *   n_lambda/n_alpha - number of hyperparameters for each list
 *   nfolds           - fold count
 *   thresh           - coordinate-descent threshold
 *   maxit            - maximum number of sweeps for coordinate descent
 *   lambda_sel/alpha_sel - pointer to the lambda/alpha value to use, or if *lambda_sel < 0 then run CV. 
 *                          After completion lambda_sel/alpha_sel is a pointer to the chosen values
 *   noise_col - length-n, added noise to final prediction for missing values for MICE DURR
 *   beta_out  - This target station's outputted coefficients  
 *   boot_col  - length-n int; bootstrapped sample drawn with replacement from Z where target is non-missing
 * 
 *                          
 * Returns 0 on success, otherwise returns:
 *   1 - for a memory allocation error
 *   2 - for invalid inputs
 *   3 - if maxit was reached by coord_descent_cov()
 *   4 - data cannot support a fit, either a constant station or too few observations per fold is breaking standardise()
 * 
 *
 */                       
static int column_impute(const double *df_old, double *target_out, int n, int p, int col, const int *missing_col, const int *fold_col,
                         const double *lambdas, int n_lambda, const double *alphas,  int n_alpha, int nfolds, double thresh, int maxit,
                         double *lambda_sel, double *alpha_sel, const double *noise_col, double *beta_out, const int *boot_col)
{
    const int p_cov = p - 1;   // number of other stations

    /* Initialise target_out by copying the target column from df_old into target_out */
    memcpy(target_out, df_old + (size_t)col * n, (size_t)n * sizeof(double));

    /* Count observed/missing rows for this station */
    int n_obs = 0;
    int n_miss = 0;
    for (int i = 0; i < n; ++i) {
        if (missing_col[i]) ++n_miss; 
        else ++n_obs;
    }

    if (n_miss == 0) { // return early if there's no missing entries
        for (int k = 0; k < p; ++k) beta_out[k] = 0.0;
        return 0; 
    }

    /* Allocating memory for all arrays */
    int *miss_days = malloc((size_t)n_miss * sizeof(int));            // missing indices of this station
    double *x_obs  = malloc((size_t)n_obs  * p_cov * sizeof(double)); // observed entries of other stations
    double *x_miss = malloc((size_t)n_miss * p_cov * sizeof(double)); // missing entries of other stations
    double *y_obs  = malloc((size_t)n_obs * sizeof(double));          // observed values of this station
    double *beta   = malloc((size_t)p_cov * sizeof(double));          // holds fitted coefficients
    double *C      = malloc((size_t)p_cov * p_cov * sizeof(double));  // covariance matrix
    char *has_row = malloc((size_t)p_cov * sizeof(char));             // flags for previously constructed rows of C

    /* Malloc error check*/
    if (!miss_days || !x_obs || !x_miss || !y_obs || !beta || !C || !has_row) {
        free(miss_days); free(x_obs); free(x_miss); free(y_obs); free(beta);
        free(C); free(has_row);
        return 1;
    }

    /* Split entries of target station into observed/missing */
    int n_im = 0;
    for (int i = 0; i < n; ++i) {
        if (missing_col[i]) miss_days[n_im++] = i;
    }

    /* Constructing each covariate station by filling values into x_obs and x_miss depending on which day the target was observed/missing */
    for (int k = 0; k < p_cov; ++k) { 

        /* Loop through all stations except the target */
        int station = (k < col) ? k : k + 1; // for columns before target station k, for columns after target station k+1, not including target index k
        
        const double *station_col = df_old + (size_t)station * n; // pointer to the start of this covariates column in df_old
        double *x_obs_col = x_obs + (size_t)k * n_obs;            // pointer to this covariate's column of x_obs
        double *x_miss_col = x_miss + (size_t)k * n_miss;         // pointer to this covariate's column of x_miss


        int constant = 1;
        /* Fill this covariate station's values on the resampled observed days into x_obs */
        for (int io = 0; io < n_obs; ++io) {
            int day = boot_col[io] - 1;         // day index of the io-th resampled day (boot_col from R)
            x_obs_col[io] = station_col[day];
            if (x_obs_col[io] != x_obs_col[0]) constant = 0; // flag if the array non-constant 
        }

        /* Error check for a constant (all the same values) covariate station across every resampled day (breaks standardise function) */
        if (constant == 1) {
            free(miss_days); free(x_obs); free(x_miss); free(y_obs); free(beta); free(C); free(has_row);
            return 4;
        }

        /* Fill this covariate station's values on the days the target station was missing into x_miss */
        for (int im = 0; im < n_miss; ++im) {
            int day = miss_days[im];            // day index of the im-th missing day
            x_miss_col[im] = station_col[day];
        }
    }

    int constant = 1;
    /* Construct the target vector y_obs from the resampled observed days */
    const double *target_col = df_old + (size_t)col * n; // pointer to the start of the target station's column in df_old
    for (int io = 0; io < n_obs; ++io) {
        int day = boot_col[io] - 1; // day index of the io-th resampled day
        y_obs[io] = target_col[day];
        if (y_obs[io] != y_obs[0]) constant = 0; // flag if the array non-constant  
    }

    /* Error check for constant (all the same values) target station across every resampled day (breaks standardise function) */
    if (constant) {
        free(miss_days); free(x_obs); free(x_miss); free(y_obs); free(beta); free(C); free(has_row);
        return 4;
    }


    double lambda = *lambda_sel;
    double alpha  = *alpha_sel;

    /* Cross-validation only if the hyperparameters were not already chosen */
    if (lambda < 0.0) {

        double *cvm   = calloc((size_t)n_alpha * n_lambda, sizeof(double)); // Mean CV error accumulated over all folds, with one entry for each lambda, alpha pair
        int *fold_obs = malloc((size_t)n_obs * sizeof(int)); // CV fold id for each observed day
        int *train_id = malloc((size_t)n_obs * sizeof(int)); // observed days used to fit this fold
        int *val_id   = malloc((size_t)n_obs * sizeof(int)); // observed days held out for validation for this fold
        /* xt, xv, yt and yv are of size n_obs, since the their sizes change for each fold */
        double *xt = malloc((size_t)n_obs * p_cov * sizeof(double)); // covariate values on the training days
        double *xv = malloc((size_t)n_obs * p_cov * sizeof(double)); // covariate values on the validation days
        double *yt = malloc((size_t)n_obs * sizeof(double)); // target values on the training days
        double *yv = malloc((size_t)n_obs * sizeof(double)); // target values on the validation days
        double *beta_path = malloc((size_t)p_cov * n_lambda * sizeof(double)); // fitted beta coefficients at every lambda
        double *int_path  = malloc((size_t)n_lambda * sizeof(double)); // fitted intercept at every lambda

        /* Malloc error check */
        if (!cvm || !fold_obs || !train_id || !val_id || !xt || !xv || !yt || !yv || !beta_path || !int_path ) {
            free(cvm); free(fold_obs); free(train_id); free(val_id); free(xt); free(xv); free(yt); free(yv);free(miss_days); free(x_obs); free(x_miss); free(y_obs); 
            free(beta); free(beta_path); free(int_path); free(C); free(has_row);
            return 1;
        }

        /* Check if there's enough data for CV and also more than 2 folds, otherwise return invalid inputs */
        if (n_obs < nfolds || nfolds < 2) {
            free(cvm); free(fold_obs); free(train_id); free(val_id); free(xt); free(xv); free(yt); free(yv); free(miss_days); free(x_obs);
            free(x_miss); free(y_obs); free(beta); free(beta_path); free(int_path); free(C); free(has_row);
            return 2;
        }

        /* Filling fold_obs with the fold ids for each value in fold_col */
        for (int io = 0; io < n_obs; ++io) fold_obs[io] = fold_col[io];

        /* Loop through folds */
        for (int f = 1; f <= nfolds; ++f) {

            /* Split into train and test */
            int n_train = 0;
            int n_val = 0;
            for (int io = 0; io < n_obs; ++io) {     // loop over observed entries
                if (fold_obs[io] == f) val_id[n_val++] = io;  // fill validation set indices
                else train_id[n_train++] = io;                // fill training set indices
            }

            /* Error check if the training or validation set is empty for some reason */
            if (n_train == 0 || n_val == 0) {
                free(cvm); free(fold_obs); free(train_id); free(val_id); free(xt); free(xv); free(yt); free(yv); free(miss_days); free(x_obs);
                free(x_miss); free(y_obs); free(beta); free(beta_path); free(int_path); free(C); free(has_row);
                return 2;
            }

            /* Copy each covariate station's observed values into the training and validation matrices */
            for (int k = 0; k < p_cov; ++k) {   // looping through each covariate station k 

                const double *x_obs_col = x_obs + (size_t)k * n_obs; // pointer to this covariate's column of x_obs
                double *xt_col = xt + (size_t)k * n_train;           // pointer to this covariate's column of xt
                double *xv_col = xv + (size_t)k * n_val;             // pointer to this covariate's column of xv

                /* Fill this covariate's values on the days used for fitting */
                for (int t = 0; t < n_train; ++t) xt_col[t] = x_obs_col[train_id[t]];
                /* Fill this covariate's values on the days held out for validation */
                for (int v = 0; v < n_val; ++v) xv_col[v] = x_obs_col[val_id[v]];
            }

            /* Split the target station's observed values across the same two sets of days */
            for (int t = 0; t < n_train; ++t) yt[t] = y_obs[train_id[t]];
            for (int v = 0; v < n_val;   ++v) yv[v] = y_obs[val_id[v]];

            /* Reset has_row for the new CV fold*/
            memset(has_row, 0, (size_t)p_cov * sizeof(char));

            /* Loop through each alpha option */
            for (int a = 0; a < n_alpha; ++a) {

                /* Fit the elastic net model for each lambda for this alpha */
                int status = elastic_net_cov(xt, yt, n_train, p_cov, lambdas, n_lambda, alphas[a], thresh, maxit, beta_path, int_path, C, has_row);

                /* Free everything if function fails and return error code */
                if (status != 0) {
                    free(cvm); free(fold_obs); free(train_id); free(val_id); free(xt); free(xv); free(yt); free(yv); free(miss_days); free(x_obs);
                    free(x_miss); free(y_obs); free(beta); free(beta_path); free(int_path); free(C); free(has_row);
                    return status;
                }


                /* For each lambda, compute the MSE on the held-out days, with each fold contributing 1/nfolds of the mean */
                for (int s = 0; s < n_lambda; ++s) {
                    const double *b = beta_path + (size_t)s * p_cov;   // pointer to the coefficients fitted at the s-th lambda option in the list
                    cvm[a + (size_t)s * n_alpha] += predict_mse(xv, n_val, p_cov, b, int_path[s], yv) / nfolds;  // cvm has one entry for each lambda/alpha option
                }
            }
        }

        /* Finding the lambda,alpha pair with the lowest MSE */
        int best_a = 0;  // best alpha
        int best_l = 0;  // best lambda
        double best = cvm[0];
        for (int s = 0; s < n_lambda; ++s) {
            for (int a = 0; a < n_alpha; ++a) {

                double v = cvm[a + (size_t)s * n_alpha]; // unpack each entry of cvm

                /* If v is Nan or infinite, return error */
                if (!isfinite(v)) {
                    free(cvm); free(fold_obs); free(train_id); free(val_id); free(xt); free(xv); free(yt); free(yv); free(miss_days); free(x_obs);
                    free(x_miss); free(y_obs); free(beta); free(beta_path); free(int_path); free(C); free(has_row);
                    return 4;
                }

                if (v < best) { best = v; best_a = a; best_l = s;}
            }
        }
        /* Selecting the best alpha and lambda from the list */
        alpha  = alphas[best_a];
        lambda = lambdas[best_l];

        free(cvm); free(fold_obs); free(train_id); free(val_id); free(xt); free(xv); free(yt); free(yv); free(beta_path); free(int_path);
    }

    /* reset has_row for the final fit */
    memset(has_row, 0, (size_t)p_cov * sizeof(char));

    /* Final fit on all observed data */
    double intercept;
    int status = elastic_net_cov(x_obs, y_obs, n_obs, p_cov, &lambda, 1, alpha, thresh, maxit, beta, &intercept, C, has_row);

    /* Free everything if function fails and return error code */
    if (status != 0) {
        free(miss_days); free(x_obs); free(x_miss); free(y_obs); free(beta); free(C); free(has_row);
        return status;
    }


    /* We need the regression error for sampling */
    double s2hat = predict_mse(x_obs, n_obs, p_cov, beta, intercept, y_obs);

    /* Predict at the missing days */
    for (int im = 0; im < n_miss; ++im) { 
        double yhat = intercept;
        for (int j = 0; j < p_cov; ++j) {
            yhat += x_miss[(size_t)j * n_miss + im] * beta[j]; // Equation 2.2.1
        }
        /* Sample imputed values using DURR */
        target_out[miss_days[im]] = yhat  + sqrt(s2hat) * noise_col[im];
    }

    /* Record the hyperparameters selected for the final fit */
    *lambda_sel = lambda;
    *alpha_sel  = alpha;

    /* Record the final beta parameters */
    beta_out[0] = intercept;
    for (int k = 0; k < p_cov; ++k) beta_out[1 + k] = beta[k];

    free(miss_days); free(x_obs); free(x_miss); free(y_obs); free(beta); free(C); free(has_row);

    return 0;
}


/* A single chained equations cycle over each target station (with synchronous updates)
 *
 * Every station is imputed from df_old alone and writes its own column of df_new
 *
 *   df_old           - n x p data matrix holding the values at the start of the cycle
 *   df_new           - n x p output matrix, gets filled with the imputed values at the end of the cycle
 *   missing_idx      - n x p, nonzero entries were originally missing
 *   folds            - n x p, CV fold id for each observed day of each station
 *   n                - number of days
 *   p                - number of stations
 *   lambdas/alphas   - CV hyperparameter options list
 *   n_lambda/n_alpha - number of hyperparameters for each list
 *   nfolds           - fold count
 *   thresh           - coordinate descent threshold
 *   maxit            - maximum number of sweeps for coordinate descent
 *   nthreads         - OpenMP thread count (0 uses avaliable number of cores)
 *   lambda_sel/alpha_sel - pointer to the lambda/alpha value to use, or if *lambda_sel < 0 then run CV. 
 *                          After completion lambda_sel/alpha_sel is a pointer to the chosen values
 *   noise     - n x p, added noise to final prediction for missing values for MICE DURR
 *   betas_out - p# x p, this target station's outputted coefficients  
 *   boot      - n x p; bootstrapped sample drawn with replacement from Z where target is non-missing
 *
 * Returns 0 on success, otherwise:
 *   1 - for a memory allocation error
 *   2 - for invalid inputs
 *   3 - if maxit was reached by coord_descent_cov()
 *   4 - data cannot support a fit, either a constant station or too few observations per fold is breaking standardise()
 *   5 - df_old contains Nan or infinity values
 *   6 - if a thread fails to report a success or an errorcode
 */
int mice_impute_cycle(const double *df_old, double *df_new, const int *missing_idx, const int *folds, int n, int p, const double *lambdas, int n_lambda,
                       const double *alphas, int n_alpha, int nfolds, double thresh, int maxit, int nthreads, double *lambda_sel, double *alpha_sel,
                       const double *noise, double *betas_out, const int *boot)
{
    /* Error check for Nan or infinity values in df_old*/ 
    for (size_t i = 0; i < (size_t)n * p; ++i) {
        if (!isfinite(df_old[i])) return 5;
    }


    if (nthreads <= 0) nthreads = omp_get_max_threads();   /* Use all cores if nthreads is not passed */

    /* Allocate space for one error code per station */
    int *status = malloc((size_t)p * sizeof(int)); 
    if (!status) return 1; // memory allocation error check
    for (int col = 0; col < p; ++col) status[col] = 6; // set default equal to 6 incase column_impute doesn't return for a target station

    /* For each target station, in parallel (synchronous updates), fit an elastic net model to predict missing values */
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) // OpenMP distributes this loop across threads, dynamic scheduling was found to be fastest
    for (int col = 0; col < p; ++col) {

        /* Call column_impute for each target station, return the error code  */
        status[col] = column_impute(df_old, df_new + (size_t)col * n, n, p, col, missing_idx + (size_t)col * n, folds + (size_t)col * n,
                      lambdas, n_lambda, alphas, n_alpha, nfolds, thresh, maxit, &lambda_sel[col], &alpha_sel[col], noise + (size_t)col*n, betas_out + (size_t)col*p, boot + (size_t)col*n);
    }

    int failed = 0;
    for (int col = 0; col < p; ++col) {
        /* If any station returns an error from column_impute, return the first error code in the list */
        if (status[col] != 0) { failed = status[col]; break; }
    }

    free(status);

    return failed; // if successful, failed should be 0 meaning success
}
