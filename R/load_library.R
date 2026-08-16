# Filename: load_library.R
#
# Description: Loads the compiled C shared library (ENCE_MICE_DURR_impute.so) 
# and exposes R functions that wrap the underlying C implementation via .Call()
#
# Compile C code in src directory to create .so file
# 
# 
# Author: Marc Lane
# Version: 6.0
# Date: 2026-08-16

# Find the shared library ENCE_MICE_DURR_impute.so in src directory
library_path <- file.path("src", paste0("ENCE_MICE_DURR_impute", .Platform$dynlib.ext))

# Error check if failed to find .so file
if (!file.exists(library_path)) {
  stop("Shared library not found at ", library_path,". Compile C code in src directory using build.sh")
}

# Load the shared library
dyn.load(library_path)

ence_impute_cycle <- function(df_old, missing_idx, folds, lambda_sel, alpha_sel,  lambda_options, alpha_options, 
                              nfolds, thresh, maxit, nthreads, cycle) {

  df_matrix      <- as.matrix(df_old)
  missing_matrix <- as.matrix(missing_idx)
  folds_matrix   <- as.matrix(folds)
  n <- nrow(df_matrix)
  p <- ncol(df_matrix)

  # Error check for valid df_old
  if (!is.numeric(df_matrix)) {
    stop("ence_impute_cycle: df_old must contain only numeric columns")
  }

  # Make sure matrices are correct dimensions
  if (nrow(missing_matrix) != n || ncol(missing_matrix) != p) {
    stop("ence_impute_cycle: missing_matrix is ", nrow(missing_matrix), " by ", ncol(missing_matrix)," but must match df_old at ", n, " by ", p)
  }
  if (nrow(folds_matrix) != n || ncol(folds_matrix) != p) {
    stop("ence_impute_cycle: folds_matrix is ", nrow(folds_matrix), " by ", ncol(folds_matrix)," but must match df_old at ", n, " by ", p)
  }

  # Make sure CV options are correct lenght
  if (length(lambda_sel) != p) {
    stop("ence_impute_cycle: lambda_sel has length ", length(lambda_sel), " but must have one entry per station, so length ", p)
  }
  if (length(alpha_sel) != p) {
    stop("ence_impute_cycle: alpha_sel has length ", length(alpha_sel), " but must have one entry per station, so length ", p)
  }

  # Make sure lambdas are given in descending order for effective warmstarting
  if (is.unsorted(rev(lambda_options))) {
    stop("ence_impute_cycle: lambda_options must be in descending order for the warm start")
  }

  storage.mode(df_matrix)      <- "double"
  storage.mode(missing_matrix) <- "integer"
  storage.mode(folds_matrix)   <- "integer"

  lam <- as.double(lambda_sel); 
  lam[is.na(lam)] <- -1.0
  alp <- as.double(alpha_sel); 
  alp[is.na(alp)] <- -1.0

  .Call("ence_impute_cycle_R", df_matrix, missing_matrix, folds_matrix, as.double(lambda_options), as.double(alpha_options), as.integer(nfolds), as.double(thresh), 
                               as.integer(maxit), as.integer(nthreads), as.integer(cycle), lam, alp, PACKAGE = "ENCE_MICE_DURR_impute")
}

mice_impute_cycle <- function(df_old, missing_idx, folds, lambda_sel, alpha_sel, boot, noise, 
                              lambda_options, alpha_options, nfolds, thresh, maxit, nthreads) {
  
  df_matrix      <- as.matrix(df_old)
  missing_matrix <- as.matrix(missing_idx)
  folds_matrix   <- as.matrix(folds)
  boot_matrix    <- as.matrix(boot)
  noise_matrix   <- as.matrix(noise)

  n <- nrow(df_matrix)
  p <- ncol(df_matrix)

  # Error check for valid df_old
  if (!is.numeric(df_matrix)) {
    stop("mice_impute_cycle: df_old must contain only numeric columns")
  }

  # Make sure matrices are correct dimensions
  if (nrow(missing_matrix) != n || ncol(missing_matrix) != p) {
    stop("mice_impute_cycle: missing_matrix is ", nrow(missing_matrix), " by ", ncol(missing_matrix)," but must match df_old at ", n, " by ", p)
  }
  if (nrow(folds_matrix) != n || ncol(folds_matrix) != p) {
    stop("mice_impute_cycle: folds_matrix is ", nrow(folds_matrix), " by ", ncol(folds_matrix)," but must match df_old at ", n, " by ", p)
  }
  if (nrow(boot_matrix) != n || ncol(boot_matrix) != p) {
    stop("mice_impute_cycle: boot_matrix is ", nrow(boot_matrix), " by ", ncol(boot_matrix)," but must match df_old at ", n, " by ", p)
  }
  if (nrow(noise_matrix) != n || ncol(noise_matrix) != p) {
    stop("mice_impute_cycle: noise_matrix is ", nrow(noise_matrix), " by ", ncol(noise_matrix)," but must match df_old at ", n, " by ", p)
  }



  # Make sure CV options are correct lenght
  if (length(lambda_sel) != p) {
    stop("mice_impute_cycle: lambda_sel has length ", length(lambda_sel), " but must have one entry per station, so length ", p)
  }
  if (length(alpha_sel) != p) {
    stop("mice_impute_cycle: alpha_sel has length ", length(alpha_sel), " but must have one entry per station, so length ", p)
  }

  # Make sure lambdas are given in descending order for effective warmstarting
  if (is.unsorted(rev(lambda_options))) {
    stop("mice_impute_cycle: lambda_options must be in descending order for the warm start")
  }

  storage.mode(df_matrix)      <- "double"
  storage.mode(missing_matrix) <- "integer"
  storage.mode(folds_matrix)   <- "integer"
  storage.mode(boot_matrix)    <- "integer"
  storage.mode(noise_matrix)   <- "double"

  lam <- as.double(lambda_sel); 
  lam[is.na(lam)] <- -1.0
  alp <- as.double(alpha_sel); 
  alp[is.na(alp)] <- -1.0

  .Call("mice_impute_cycle_R", df_matrix, missing_matrix, folds_matrix, as.double(lambda_options), as.double(alpha_options),
        as.integer(nfolds), as.double(thresh), as.integer(maxit),as.integer(nthreads), lam, alp, noise_matrix, boot_matrix, PACKAGE = "ENCE_MICE_DURR_impute")
}