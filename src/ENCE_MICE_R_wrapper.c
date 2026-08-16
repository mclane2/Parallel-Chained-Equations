
/**
 * Filename: ENCE_MICE_R_wrapper.c
 * 
 * Description: This is the R SEXP wrapper so that R can call ence_impute_cycle and mice_impute_cycle implemented in C.
 * Unpacks R objects into plain C types, runs the C functions, and packs results back into an R list.
 * 
 * Author: Marc Lane
 * Version: 5.0
 * Date: 2026-08-15
 */
#include <R.h>
#include <Rinternals.h>
#include "ENCE_impute.h"
#include "MICE_DURR_impute.h"

/* ENCE */
SEXP ence_impute_cycle_R(SEXP df_old_sexp, SEXP missing_sexp, SEXP folds_sexp, SEXP lambdas_sexp, SEXP alphas_sexp, SEXP nfolds_sexp, 
                         SEXP thresh_sexp, SEXP maxit_sexp, SEXP nthreads_sexp,  SEXP cycle_sexp, SEXP lambda_sel_sexp, SEXP alpha_sel_sexp)
{
    int n = nrows(df_old_sexp);
    int p = ncols(df_old_sexp);

    /* Output matrix df_new */
    SEXP df_new_sexp = PROTECT(allocMatrix(REALSXP, n, p));

    /* Duplicate the in/out hyperparameter vectors */
    SEXP lambda_out = PROTECT(duplicate(lambda_sel_sexp));
    SEXP alpha_out  = PROTECT(duplicate(alpha_sel_sexp));

    int failed = ence_impute_cycle(REAL(df_old_sexp), REAL(df_new_sexp), INTEGER(missing_sexp), INTEGER(folds_sexp), n, p, REAL(lambdas_sexp), LENGTH(lambdas_sexp),
                                   REAL(alphas_sexp), LENGTH(alphas_sexp), asInteger(nfolds_sexp), asReal(thresh_sexp), asInteger(maxit_sexp), asInteger(nthreads_sexp),
                                   asInteger(cycle_sexp), REAL(lambda_out), REAL(alpha_out));
    /* Translate the C error code into an R error message */
    if (failed == 1) Rf_error("ence_impute_cycle: memory allocation failed");
    if (failed == 2) Rf_error("ence_impute_cycle: invalid inputs passed");
    if (failed == 3) Rf_error("ence_impute_cycle: coordinate descent reached maxit without converging");
    if (failed == 4) Rf_error("ence_impute_cycle: data cannot support a fit (constant station, or too few observations per fold breaks standardise() function)");
    if (failed == 5) Rf_error("ence_impute_cycle: df_old contains NaN or infinite values (either present in the input or produced by the previous cycle's fit)");
    if (failed == 6) Rf_error("ence_impute_cycle: a thread failed to report a status code");
    if (failed == 7) Rf_error("ence_impute_cycle: could not open betas.csv for writing");
    if (failed)      Rf_error("ence_impute_cycle: unrecognised error code %d", failed);

    /* Pack into a named list: df, lambda, alpha  */
    SEXP result = PROTECT(allocVector(VECSXP, 3));
    SEXP names  = PROTECT(allocVector(STRSXP, 3));

    /* Fill list */
    SET_VECTOR_ELT(result, 0, df_new_sexp);
    SET_VECTOR_ELT(result, 1, lambda_out);
    SET_VECTOR_ELT(result, 2, alpha_out);
    /* Attach labels */
    SET_STRING_ELT(names, 0, mkChar("df"));
    SET_STRING_ELT(names, 1, mkChar("lambda"));
    SET_STRING_ELT(names, 2, mkChar("alpha"));
    /* Attach label to the list */
    setAttrib(result, R_NamesSymbol, names);

    // Unprotect the allocated vectors
    UNPROTECT(5);
    return result;
}


/* MICE DURR */
SEXP mice_impute_cycle_R(SEXP df_old_sexp, SEXP missing_sexp, SEXP folds_sexp, SEXP lambdas_sexp, SEXP alphas_sexp,
                         SEXP nfolds_sexp, SEXP thresh_sexp, SEXP maxit_sexp, SEXP nthreads_sexp, SEXP lambda_sel_sexp, SEXP alpha_sel_sexp,
                         SEXP noise_sexp, SEXP boot_sexp)
{
    int n = nrows(df_old_sexp);
    int p = ncols(df_old_sexp);

    /* Output matrix df_new */
    SEXP df_new_sexp = PROTECT(allocMatrix(REALSXP, n, p));

    /* Output betas */
    SEXP betas_sexp = PROTECT(allocMatrix(REALSXP, p, p));

    /* Duplicate the in/out hyperparameter vectors */
    SEXP lambda_out = PROTECT(duplicate(lambda_sel_sexp));
    SEXP alpha_out  = PROTECT(duplicate(alpha_sel_sexp));

    int failed = mice_impute_cycle(REAL(df_old_sexp), REAL(df_new_sexp), INTEGER(missing_sexp), INTEGER(folds_sexp), n, p, REAL(lambdas_sexp), LENGTH(lambdas_sexp),
                                   REAL(alphas_sexp), LENGTH(alphas_sexp), asInteger(nfolds_sexp), asReal(thresh_sexp), asInteger(maxit_sexp), asInteger(nthreads_sexp),
                                   REAL(lambda_out), REAL(alpha_out), REAL(noise_sexp), REAL(betas_sexp), INTEGER(boot_sexp));
    /* Translate the C error code into an R error message */
    if (failed == 1) Rf_error("mice_impute_cycle: memory allocation failed");
    if (failed == 2) Rf_error("mice_impute_cycle: invalid inputs passed");
    if (failed == 3) Rf_error("mice_impute_cycle: coordinate descent reached maxit without converging");
    if (failed == 4) Rf_error("mice_impute_cycle: data cannot support a fit (constant station, or too few observations per fold breaks standardise() function)");
    if (failed == 5) Rf_error("mice_impute_cycle: df_old contains NaN or infinite values");
    if (failed == 6) Rf_error("mice_impute_cycle: a thread failed to report a status code");
    if (failed)      Rf_error("mice_impute_cycle: unrecognised error code %d", failed);

    /* Pack into a named list: df, lambda, alpha, betas */
    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SEXP names  = PROTECT(allocVector(STRSXP, 4));

    /* Fill list */
    SET_VECTOR_ELT(result, 0, df_new_sexp);
    SET_VECTOR_ELT(result, 1, lambda_out);
    SET_VECTOR_ELT(result, 2, alpha_out);
    SET_VECTOR_ELT(result, 3, betas_sexp);

    /* Attach labels */
    SET_STRING_ELT(names, 0, mkChar("df"));
    SET_STRING_ELT(names, 1, mkChar("lambda"));
    SET_STRING_ELT(names, 2, mkChar("alpha"));
    SET_STRING_ELT(names, 3, mkChar("betas"));
    /* Attach label to the list */
    setAttrib(result, R_NamesSymbol, names);

    // Unprotect the allocated vectors
    UNPROTECT(6);
    return result;
}