# Parallel Chained Equations

This repository showcases the work from my M.Sc. thesis project, *Chained Equations for Incomplete Climate Networks*, supervised by Brian O'Sullivan (Met Éireann), for completion of my M.Sc. in High Performance Computing at Trinity College Dublin.

The problem of missing data is a pervasive feature of many areas of empirical research. Imputation is the infilling of missing values with plausible estimates that preserve the statistical properties of the underlying data-generating process, rather than simply producing a complete-looking dataset. Imputation is particularly important for climate monitoring networks due to their high levels of missingness, since continuous long-term climate records are essential for producing forecasts and analysing long-term climate change.

([O'Sullivan and Kelly, 2024](https://doi.org/10.1002/joc.8513)) introduced two new imputation algorithms, elastic net chained equations (ENCE) and multiple imputation by chained equations with direct use of regularized regression (MICE DURR). Both methods were originally developed for imputing monthly rainfall data from Met Éireann's network throughout the Republic of Ireland. The algorithms use elastic net models to deal with the high dimensionality of Met Éireann's monthly rainfall network, where the number of stations exceeded the number of observations. These methods were constructed using the chained equations framework, which fits a model for each station in sequence, cycling through the entire dataset until the imputed values stop changing.

My thesis focused on imputing daily rainfall data from Met Éireann's monitoring network. Moving from monthly to daily resolution makes both methods far more expensive to run, which turns imputation at this scale into a high performance computing problem. To improve the performance of ENCE and MICE DURR, my thesis introduced synchronous updates, which is a modification to the structure of chained equations that removes a serial dependency and lets the models within a cycle be fitted in parallel. Synchronous updates can be implemented in R using the parallel package with minimal changes to the original code, leading to a substantial speedup. I also implemented both methods in C with OpenMP for shared memory parallelism. Moving the chained equations section of each algorithm to C required writing the cyclical coordinate descent solver for the elastic net models in C. This allowed me to add SIMD vectorisation pragmas on the most computationally expensive sections of the solver to obtain a further speedup. As well as parallelising synchronous updates across cores, MICE DURR uses multiple imputation which can be trivially parallelised across the nodes of a cluster. This repository presents MICE DURR with an OpenMP-MPI hybrid structure which uses synchronous updates parallelised across cores and multiple imputation parallelised across nodes. My thesis also compared the performance of ENCE and MICE DURR to inverse distance weighting (IDW), multiple linear regression (MLR), MissForest, probabilistic principal component analysis (PPCA) and scikit-learn's IterativeImputer. The C implementations of ENCE and MICE DURR showcased in this repository were found to be the most accurate of all methods tested. 

This repository presents three implementations of both ENCE and MICE DURR: the original serial versions written entirely in R, the synchronous versions also in R, and the C implementation with synchronous updates parallelised across threads using OpenMP. This repository uses simulated daily data (data/simulated_daily_rainfall.csv), which simulates 1 year of daily rainfall in Ireland across 200 stations. The user can edit data/simulate_rainfall.Rmd to create a different dataset. The real Met Éireann daily rainfall dataset is private and requires permission to access from Met Éireann. 

The original serial implementations of ENCE and MICE DURR, and the data simulation file were adopted from Brian O'Sullivan's repository [ENCE_MICE_DURR](https://github.com/BrianOSullivan-2000/ENCE_MICE_DURR). 

To access the working directory with the code used to obtain the results for my M.Sc. thesis, see [HPC-Thesis](https://github.com/mclane2/HPC-Thesis).

## Results

The table below shows speedups compared to the original, serial R implementation from demo.Rmd on the simulated dataset included in this repository (200 stations, 365 days). The figures here were obtained on my laptop with 6 cores:

| Method | Implementation| Speedup vs serial | Speedup vs R synchronous |
|---|---|---|---|---|
| ENCE | R synchronous | 3.1 x | — |
| ENCE | C / OpenMP | 11.4 x | 3.6 x |
| MICE DURR (m=1) | R synchronous | 2.6 x | — |
| MICE DURR (m=1)| C / OpenMP |  7 x | 2.7 x |

The figures here differ from those presented in my thesis, which was for imputing ten years of real daily rainfall data across 523 stations. The table below shows speedups from my thesis. My thesis figures were obtained on the Seagull cluster operated by TCHPC at Trinity College Dublin. Each node on Seagull has 16 cores (two Intel Xeon E5-2630 v3 CPUs at 2.40GHz) with 64GB of RAM:

| Method | Implementation | Speedup vs serial | Speedup vs R synchronous |
|---|---|---|---|---|
| ENCE | R synchronous |6.46 x | — |
| ENCE | C / OpenMP | 8.62 x | 1.33 x |
| MICE DURR (m=4) | R synchronous | 4.80 x | — |
| MICE DURR (m=4) | C / OpenMP |5.71 x | 1.19 x |

The ENCE and MICE DURR took much longer for the results presented in my thesis, which was over ten years compared to a single year for the simulated dataset. The larger speedup of the C implementation for the simulated dataset can be explained by R's overhead getting amortised for larger datasets. The benefit of synchronous updates is clearest on the ten year dataset, where they account for most of the total speedup in both methods.

## Requirements

**To run the demo**

**To re-simulate the data**

## Running the demo

## Repository layout

## Data

## Citation
