# Parallel Chained Equations

This repository showcases the work from my thesis project, *Chained Equations for Incomplete Climate Networks*, supervised by Brian O'Sullivan (Met Éireann), and completed as part of my M.Sc. in High Performance Computing at Trinity College Dublin.

Missing data is a pervasive problem across many areas of research. Imputation is the infilling of missing values with plausible estimates that preserve the statistical properties of the underlying data-generating process, rather than simply producing a complete-looking dataset. Imputation is particularly important for climate monitoring networks due to their high levels of missingness, since continuous long-term climate records are essential for producing forecasts and analysing long-term climate change.

[O'Sullivan and Kelly, 2024](https://doi.org/10.1002/joc.8513) introduced two new imputation algorithms, elastic net chained equations (ENCE) and multiple imputation by chained equations with direct use of regularized regression (MICE DURR). Both methods were originally developed for imputing monthly rainfall data from Met Éireann's network throughout the Republic of Ireland. The algorithms use elastic net models to deal with the high dimensionality of Met Éireann's monthly rainfall network, where the number of stations exceeded the number of observations. These methods were constructed using the chained equations framework, which fits a model for each station in sequence, cycling through the entire dataset until the imputed values stop changing.

My thesis focused on imputing daily rainfall data from Met Éireann's monitoring network. Moving from monthly to daily resolution makes both methods far more expensive to run, which turns imputation at this scale into a high performance computing problem. To improve the performance of ENCE and MICE DURR, my thesis introduced synchronous updates, a modification to the structure of chained equations that removes a serial dependency and lets the models within a cycle be fitted in parallel. Synchronous updates can be implemented in R using the parallel package with minimal changes to the original code, leading to a substantial speedup. I also implemented both methods in C with OpenMP for shared memory parallelism. Moving the chained equations section of each algorithm to C required writing the cyclical coordinate descent solver for the elastic net models in C. This allowed me to add SIMD vectorisation pragmas on the most computationally expensive sections of the solver to obtain a further speedup. As well as parallelising synchronous updates across cores, MICE DURR uses multiple imputation which can be trivially parallelised across the nodes of a cluster. This repository presents MICE DURR with a doubly parallel structure which uses synchronous updates parallelised across cores and multiple imputation parallelised across nodes. My thesis also compared the performance of ENCE and MICE DURR to inverse distance weighting (IDW), multiple linear regression (MLR), MissForest, probabilistic principal component analysis (PPCA) and scikit-learn's `IterativeImputer`. The C implementations of ENCE and MICE DURR showcased in this repository were the most accurate of all methods tested. 

This repository presents three implementations of both ENCE and MICE DURR: the original serial versions written entirely in R, the synchronous versions also in R, and the C implementation with synchronous updates parallelised across threads using OpenMP. This repository uses simulated daily data (data/simulated_daily_rainfall.csv), which simulates one year of daily rainfall in Ireland across 200 stations. The user can edit data/simulate_rainfall.Rmd to create a different dataset. The real daily rainfall dataset used in my thesis is private and requires permission to access from Met Éireann. 

The original serial implementations of ENCE and MICE DURR, and the data simulation file were adapted from Brian O'Sullivan's repository [ENCE_MICE_DURR](https://github.com/BrianOSullivan-2000/ENCE_MICE_DURR). 

To access the working directory with the code used to obtain the results for my M.Sc. thesis, see [HPC-Thesis](https://github.com/mclane2/HPC-Thesis).

## Results

The table below shows speedups compared to the original, serial R implementation from `demo.Rmd` on the simulated dataset included in this repository (200 stations, 365 days). The figures here were obtained on my laptop with 6 cores:

| Method | Implementation| Speedup vs serial | Speedup vs R synchronous |
|---|---|---|---|
| ENCE | R synchronous | 5.1× | — |
| ENCE | C / OpenMP | 17.3× | 3.4× |
| MICE DURR (m=1) | R synchronous | 3.9× | — |
| MICE DURR (m=1)| C / OpenMP |  8.1× | 2.1× |

The figures here differ from those presented in my thesis, which was for imputing ten years of real daily rainfall data across 523 stations. The table below shows speedups from my thesis. My thesis figures were obtained on the Seagull cluster operated by TCHPC at Trinity College Dublin. Each node on Seagull has 16 cores (two Intel Xeon E5-2630 v3 CPUs at 2.40GHz) with 64GB of RAM:

| Method | Implementation | Speedup vs serial | Speedup vs R synchronous |
|---|---|---|---|
| ENCE | R synchronous |6.5× | — |
| ENCE | C / OpenMP | 8.6× | 1.3× |
| MICE DURR (m=4) | R synchronous | 4.8× | — |
| MICE DURR (m=4) | C / OpenMP |5.7× | 1.2× |

ENCE and MICE DURR took much longer for the results presented in my thesis because the dataset covered ten years compared to a single year for the simulated dataset. The larger speedup of the C implementation on the simulated dataset is explained by R's overhead being amortised. On the small simulated dataset that overhead is a large fraction of the runtime, so avoiding it in C gives a large speedup, while on the ten year dataset it is negligible. The benefit of synchronous updates is clearest on the ten year dataset, where they account for most of the total speedup in both methods.

## Requirements

This code was developed and tested on Ubuntu running on WSL (Windows Subsystem for Linux). It should work on any Linux system with a C compiler that supports OpenMP, but it has not been tested on macOS or native Windows.

**To run the demo**

R (≥ 4.0), install with: 
```bash
sudo apt install r-base
```
A C compiler with OpenMP support (gcc), install with:
```bash
sudo apt install build-essential
```

C code was compiled into a shared library that can be called by R using `R CMD SHLIB`, which requires: 
```bash
sudo apt install r-base-dev
```

The following R packages are needed:

| Package | Used for |
|---|---|
| `dplyr`, `tidyr` | reshaping between long and wide tables |
| `glmnet`, `glmnetUtils` | elastic net fits and cross validation in the R implementations |
| `sp`, `gstat` | inverse distance weighting for the starting values (`R/idw.R`) |
| `spacetime` | only needed if you call `IDW_ST` with `C = NULL`, which estimates the anisotropy. The demo passes `C = 15000`, so it is not used there |
| `parallel` | forked workers for synchronous updates. Should already be installed with R, no install needed |

```r
install.packages(c("dplyr", "tidyr", "glmnet", "glmnetUtils", "sp", "gstat", "spacetime"))
```

Synchronous updates in the R implementations use `parallel::mclapply`, which relies on forking processes and does not run in parallel on native Windows.

To knit `demo.Rmd`, you also need the `rmarkdown` and `knitr` R packages:
```r
install.packages(c("rmarkdown", "knitr"))
```

This may require you to install pandoc, which can be installed with:
```bash
sudo apt install pandoc
```

**To re-simulate the data**

`data/simulate_rainfall.Rmd` needs:

| Package | Used for |
|---|---|
| `sf` | reading the Ireland boundary and plotting |
| `spatstat.geom`, `spatstat.random` | sampling random station locations |
| `ggplot2` | plotting the simulated stations |

```r
install.packages(c("sf", "spatstat.geom", "spatstat.random", "ggplot2"))
```

The `sf` package may require you to install:
```bash
sudo apt install libgdal-dev libgeos-dev libproj-dev libudunits2-dev
```


## Running the demo

Clone the repository and move to it:
```bash
git clone https://github.com/mclane2/Parallel-Chained-Equations.git
cd Parallel-Chained-Equations
```

Run the demo by knitting `demo.Rmd`:
```bash
Rscript -e "rmarkdown::render('demo.Rmd')"
```
The results of `demo.Rmd` are written to `demo.html`, which can be opened on WSL using:
```bash
explorer.exe demo.html
```

Running `demo.Rmd` automatically compiles the C code in the src directory into a shared library `src/ENCE_MICE_DURR_impute.so`. To compile the C code ahead of time, you can compile it manually using 
```bash
cd src && chmod +x compile.sh  && ./compile.sh && cd ..
```

To determine the number of cores on the user's machine for synchronous updates, I used the function `parallel::detectCores(logical = FALSE)`. On WSL/Linux, this function doesn't work reliably, and will probably return the number of logical cores (my laptop has 6 cores but 12 logical cores). If you know the number of cores on your machine, I recommend manually setting `n_cores` in `demo.Rmd`.

For the simulated data (200 stations, 365 days) run on 6 cores on my laptop, `demo.Rmd` altogether takes approximately 7.5 minutes to run.

## Repository layout

```
Parallel-Chained-Equations/
├── demo.Rmd                        # Runs and compares all six implementations
├── README.md
├── R/
│   ├── ENCE_serial.R               # Original serial ENCE
│   ├── ENCE_synchronous.R          # ENCE with synchronous updates (uses R parallelism)
│   ├── ENCE_C.R                    # ENCE calling the C solver
│   ├── MICE_DURR_serial.R          # Original serial MICE DURR
│   ├── MICE_DURR_synchronous.R     # MICE DURR with synchronous updates (uses R parallelism)
│   ├── MICE_DURR_C.R               # MICE DURR calling the C solver
│   ├── idw.R                       # Inverse distance weighting for initial values
│   └── load_library.R              # Loads the shared library
├── src/
│   ├── elastic_net_functions.c     # Main cyclical coordinate descent solver functions
│   ├── elastic_net.c               # Driver for functions in elastic_net_functions.c
│   ├── ENCE_impute.c               # One synchronous ENCE chained equations cycle
│   ├── MICE_DURR_impute.c          # One synchronous MICE DURR chained equations cycle
│   ├── ENCE_MICE_R_wrapper.c       # .Call entry points to make C functions visible to R
│   ├── .h files                    # For linking all C files
│   └── compile.sh                  # Run to compile everything in src into ENCE_MICE_DURR_impute.so
└── data/
    ├── simulated_daily_rainfall.csv   # 200 stations, 365 days
    ├── counties.json                  # Ireland shapefile
    └── simulate_rainfall.Rmd          # Script used to generate simulated data
```

## Citation

Please cite the following paper if you use this code: O'Sullivan, B., & Kelly, G. (2024). Infilling of high-dimensional rainfall networks through multiple imputation by chained equations. International Journal of Climatology, 44(9), 3075–3091. https://doi.org/10.1002/joc.8513

Please contact me or Brian O'Sullivan (Met Éireann) if you want to use this code. Please reach out if you have any issues with the code.

## Contact

Email: **Marc Lane** – marc.lane.3002@gmail.com

