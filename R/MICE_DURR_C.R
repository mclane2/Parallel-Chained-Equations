
# Filename: MICE_C.R 
# Author: Brian O'Sullivan (modified by Marc Lane)

# Please refer to https://doi.org/10.1002/joc.8513 for further details 
# where the method is implemented on the Irish rainfall network

## Primary imputation function

## NOTE: The dataframe must follow a "long-table" format. i.e. each row is a 
## measurement with a spatial ID and time ID. All possible ID combinations
## must be included, missing or non-missing, and duplicate 
## combinations are not allowed

## Modifications:
## - Changed MICE_get_params function to do synchronous updates.
## - Parallelism is across synchronous updates only
## - Added RMSE tracking
## - Calls C implementation of glmnet

## lambda_options must be passed in descending order. This is what allows the warmstart optimization to work effectively in the C coordinate descent solver


MICE_DURR_C <- function(df, response = "y", m = 4, 
                      max_impute_cycles = 5, 
                      max_EM_cycles = 5, 
                      hyp_cycles = 2,
                      impute_tol = 0,
                      init_method = c("mean", "idw"),
                      spatial_id = "stno", time_id = "t",
                      EM_tol = 1,
                      nclusters = 1,
                      nthreads = 0L,
                      nfolds = 10L,
                      lambda_options = c(0.2, 0.15, 0.1, 0.05),
                      alpha_options  = c(0.1, 0.35, 0.65, 0.9),
                      thresh = 1e-7,
                      maxit = 100000L,
                      outfile = "MICE_DURR_log.txt",
                      C = 15000, x_name = "east", y_name = "north",
                      transformation = {function(x) x},
                      reverse_transformation = {function(x) x},
                      ...){
  
  # df                      -Data frame of response variable and covariates
  # max_impute_cycles       -Maximum number of cycles when getting parameters
  # max_EM_cycles           -Maximum number of cycles when getting final values
  # hyp_cycles              -Number of cycles to fit hyperparameters
  # impute_tol              -Tolerance for convergence between cycles when getting parameters
  # init_method             -Method for getting starting imputed values
  # spatial_id, time_id     -Column names for locations and times
  # EM_tol                  -Tolerance for convergence of final values 
  # nclusters               -Number of clusters (nodes) if doing parallel multiple imputation
  # nthreads                -Number of threads per node parallel updates
  # nfolds                  -Number of folds for CV
  # lambda_options          -Lambda choices for CV to pick from, must be in descending order
  # alpha_options           -Alpha choices for CV to pick from
  # thresh                  -Convergence threshold for the C coordinate descent solver
  # maxit                   -Maximum number of sweeps for the C coordinate descent solver
  # outfile                 -Location to print outputs from clusters
  
  
  # Load in necessary packages
  require(dplyr); require(tidyr)
  
  init_method <- match.arg(init_method)
  
  # Label missing entries
  df <- df %>% mutate(missing = is.na(.data[[response]]))
  original_spatial_order <- unique(df[[spatial_id]])
  original_times <- df[[time_id]]
  
  # Sort the spatial locations from least missing values to most
  # (This isn't necessary, but can improve imputation accuracy)
  df <- df %>%
    group_by(across(all_of(spatial_id))) %>%
    mutate(.n_missing = sum(missing)) %>%
    arrange(.n_missing) %>%
    dplyr::select(-.n_missing)
  
  # Get starting imputed values
  # Mean imputation (mean of each spatial location)
  if(init_method == "mean"){
    df <- df %>% group_by(across(all_of(spatial_id))) %>%
      mutate(across(all_of(response), ~replace_na(., mean(., na.rm=TRUE))))
  }
  # Impute using Inverse Distance Weighting (idw)
  # refer to functions in "idw.R"
  else if(init_method == "idw"){
    df[df$missing, ] <- IDW_ST(df[!df$missing, ], df[df$missing, ], 
                               response,  C = C, 
                               x_name = x_name, y_name = y_name, ...)
  }
  else{stop("Error: invalid starting imputation method selected")}
  
  # Transformation of response
  df[response] <- transformation(df[response])
  
  # Get wide table of missing entry locations
  missing_idx <- df %>% dplyr::select(all_of(c(time_id, spatial_id)), missing) %>%
    pivot_wider(names_from = all_of(spatial_id), 
                values_from = missing) %>%
    dplyr::select(-all_of(time_id)) %>% as.matrix()
  
  # Convert dataframe to wide table
  df <- df %>% dplyr::select(all_of(c(time_id, spatial_id, response))) %>%
    pivot_wider(names_from = all_of(spatial_id), 
                values_from = all_of(response)) %>%
    dplyr::select(-all_of(time_id))
  
  # Add character to beginning of IDs in case they are just numbers
  colnames(missing_idx) <- paste0(spatial_id, colnames(missing_idx))
  names(df) <- paste0(spatial_id, names(df))

  # Setting up Multiple Imputation
  # Replicate wide table into m copies and attach a seed to each
  df_copies <- mapply(list, 222:(222+(m-1)), rep(list(df), m), 
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  if(nclusters > 1){
    require(parallel)

    # Make clusters for parallel computing
    cl <- makeCluster(nclusters, outfile = outfile)

    # Import necessary functions and the shared library to each cluster
    clusterExport(cl, c("rmse", "mice_impute_cycle", "library_path"))
    clusterEvalQ(cl, dyn.load(library_path))

    # Get parameters for every imputation model
    results <- parLapply(cl = cl, df_copies, MICE_get_params_C,
                      missing_idx = missing_idx, max_cycles = max_impute_cycles,
                      hyp_cycles = hyp_cycles, lambda_options = lambda_options,
                      alpha_options = alpha_options, nfolds = nfolds, thresh = thresh,
                      maxit = maxit, impute_tol = impute_tol, nthreads = nthreads, ...)
    stopCluster(cl)
  }
  else{
    # This provides the option to not do parallel computing
    results <- lapply(df_copies, MICE_get_params_C,
                     missing_idx = missing_idx, max_cycles = max_impute_cycles,
                     hyp_cycles = hyp_cycles, lambda_options = lambda_options,
                     alpha_options = alpha_options, nfolds = nfolds, thresh = thresh,
                     maxit = maxit, impute_tol = impute_tol, nthreads = nthreads, ...)
  }
  # Pool the regression parameters together from each replicated data set
  params <- Reduce(`+`, lapply(results, `[[`, "params")) / m
  colnames(params) <- names(df)
  rmse_histories <- lapply(results, `[[`, "rmse_history")
  lambdas_list      <- lapply(results, `[[`, "lambdas")
  alphas_list       <- lapply(results, `[[`, "alphas")
  beta_history_list <- lapply(results, `[[`, "beta_history")
  
  # Iterator and convergence tracking
  i <- 1; old_rmse <- .Machine$double.xmax; new_rmse <- .Machine$double.xmax
  df <- as.data.frame(df)
  
  # Update missing values over multiple cycles with fixed regression parameters
  while((i <= max_EM_cycles) & ((old_rmse/new_rmse > EM_tol) | (i == 1))){
  
    # Update imputed values and convergence tracking
    past_values <- df[missing_idx]
    old_rmse <- new_rmse

    # Update imputed values
    df <- MICE_impute_fixed(df, missing_idx, params, ...)
    
    # New RMSE for convergence
    new_rmse <- rmse(past_values, df[missing_idx])
    
    print(paste0("Cycle ", i, " completed"))
    i <- i + 1
  }
  if(i == max_EM_cycles){
    print("Maximum number of cycles reached")
  }
  else{
    print(paste0("Converged at cycle: ", i))
  }
  
  # Return df to long format
  names(df) <- substr(names(df), nchar(spatial_id)+1, nchar(names(df)))
  df <- pivot_longer(df, cols = everything(), 
                     names_to = spatial_id, values_to = response)
  df <- df %>%
    mutate(across(all_of(spatial_id), 
                  ~factor(.x, levels = original_spatial_order))) %>%
    arrange(across(all_of(spatial_id)))
  
  # Add times back and undo transformation
  df[time_id] <- original_times
  df[response] <- reverse_transformation(df[response])

  df <- as.data.frame(df)
  attr(df, "rmse_histories")    <- rmse_histories
  attr(df, "lambdas_list")      <- lambdas_list
  attr(df, "alphas_list")       <- alphas_list
  attr(df, "beta_history_list") <- beta_history_list
  return(df)
}

# For each replicated dataset, calculate the regression parameters
# for each imputation model
MICE_get_params_C <- function(df, missing_idx, max_cycles, hyp_cycles, lambda_options, alpha_options, nfolds, thresh, maxit, impute_tol, nthreads, ...){
  
  # df              -Input data
  # missing_idx     -Index of all originally missing values
  # max_cycles      -Max number of cycles to update values
  # hyp_cycles      -Number of cycles to fit hyperparameters
  # lambda_options  -Lambda values to consider for elastic-net
  # alpha_options   -Alpha values to consider for elastic-net
  # impute_tol      -Tolerance for convergence between cycles
  # nfolds          -CV folds
  # nthreads        -Number of threads for synchronous parallel updates

  # lambda_options MUST be passed in descending order
  # (otherwise it will break warmstart C solver)

  # Extract the random seed and data
  seed <- df[[1]]
  df <- df[[2]]

  # Set random seed
  set.seed(seed)

  # Track RMSE per cycle
  rmse_history <- c()
  beta_history <- list()

  # If there's only one option, don't do CV, use fixed alpha/lambda
  # Otherwise NA_real_ will trigger CV of lambda/allpha_options
  L0 <- if (length(lambda_options) == 1) lambda_options else NA_real_
  A0 <- if (length(alpha_options)  == 1) alpha_options  else NA_real_
  # ls/as will hold the selected lambda/alpha for each station
  ls <- rep(L0, ncol(df))
  as <- rep(A0, ncol(df))

  # Prepare bootstrapped row-index/noise/fold matrices for C implementation
  col_names <- names(df)
  n <- nrow(df); p <- ncol(df)
  boot  <- matrix(0L, n, p)    # row indices of bootstrap resample
  noise <- matrix(0,  n, p)    # Precomputed N(0,1) draws for DURR noise
  folds <- matrix(0L, n, p)    # CV fold id per bootstrap position
  # Loop through each column to construct boot and noise
  for (j in seq_len(p)) {

    ry  <- missing_idx[, j]
    obs <- which(!ry)
    # if there are no observed rows
    if (length(obs) == 0) next

    # Constructing bootstrapped columns for each station
    set.seed(seed)
    # Sample from 1,2, ... , length(obs) with replacement
    s <- sample(length(obs), length(obs), replace = TRUE)
    boot[seq_along(s), j] <- obs[s]

    # Precomputing rnorm for each column
    set.seed(seed)
    noise[seq_len(sum(ry)), j] <- rnorm(sum(ry))

    # Randomly associating a fold id to each observation
    set.seed(seed)
    folds[seq_along(s), j] <- sample(rep(1:nfolds, length.out = length(s)))
  }

  # Iterator and convergence tracking
  i <- 1; old_rmse <- .Machine$double.xmax; new_rmse <- .Machine$double.xmax
  
  # Update the imputed values over multiple cycles
  while((i <= max_cycles) & ((old_rmse/new_rmse > impute_tol) | (i == 1))){
    
    # Update imputed values and convergence tracking
    past_values <- df[missing_idx]
    old_rmse <- new_rmse
    
    # Reset lambdas and alphas if still updating hyperparameters
    if (i < hyp_cycles){
      ls <- rep(L0, ncol(df)); as <- rep(A0, ncol(df))
    }

    # Call C implementation and compute imputed values
    imputed_df <- mice_impute_cycle(df, missing_idx, folds, ls, as, boot, noise, lambda_options = lambda_options, alpha_options = alpha_options,
                                    nfolds = nfolds, thresh = thresh, maxit = maxit, nthreads = nthreads)

    beta_history[[i]] <- imputed_df$betas
    df        <- as.data.frame(imputed_df$df)
    names(df) <- col_names
    ls        <- imputed_df$lambda
    as        <- imputed_df$alpha

    new_rmse <- rmse(past_values, df[missing_idx])
    rmse_history <- c(rmse_history, new_rmse)

    print(paste0("Cycle ", i, " completed"))
    i <- i + 1
  }
  if(i == max_cycles){
    print("Maximum number of cycles reached")
  }
  else{
    print(paste0("Converged at cycle: ", i - 1))
  }

  # Final non-bootstrapped fit and all converged data with glmnet
  require(glmnet)
  betas <- matrix(0, p, p)
  
  # Get regression parameters for each column
  for (column in seq_len(p)) {

    # Skip if there's no missing data
    if (sum(missing_idx[, column]) == 0) next

    # Input matrices for glmnet
    x <- as.matrix(df[, -column])
    y <- df[[column]]

    # Fit model using elastic-net and save regression parameters
    model <- glmnet(x = x, y = y, lambda = ls[column], alpha = as[column])
    betas[, column] <- as.matrix(coef(model))[, 1]
  }

  return(list(params = betas, rmse_history = rmse_history, lambdas = ls, alphas = as, beta_history = beta_history))

}



# Update missing values with fixed regression parameters

MICE_impute_fixed <- function(df, missing_idx, params, ...){
  
  # Impute each column
  for (column in 1:ncol(df)){
    
    # Get covariates (other columns) and positions of missing values
    covariates <- as.matrix(df[, -column])
    missing_values <- missing_idx[, column]
    
    # Only care if there are missing values to be updated in column
    if(sum(missing_values) > 0){
        
      # Get fixed regression model parameters for current column
      coefs <- params[, column]
        
      # Linear model -> Imputed values calculated with X %*% coefs
      if(sum(missing_values) == 1){
        # Edge case if only one missing value in column
        X <- c(1, covariates[missing_values, ])
        df[missing_values, column] <- as.numeric(X) %*% coefs 
      }
      else{
        df[missing_values, column] <- 
          cbind(rep(1, nrow(df[missing_values, ])),
                as.matrix(covariates)[missing_values, ]) %*% coefs
      }
    }
  }
  
  return(df)
}


# Quick rmse function

rmse <- function(y1, y2){
  return(sqrt(mean((y1 - y2)^2)))
}