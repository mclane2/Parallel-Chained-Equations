#!/usr/bin/env bash
#
# Filename: compile.sh
#
# Description: Compiles the C source files into .so files (an R-loadable shared library) via R CMD SHLIB, with OpenMP
#
# To Compile input to terminal:
# chmod +x compile.sh
# ./compile.sh
#
# Author: M. Lane
# Version: 7.0 (Cleaned up for new repository)
# Date: 2026-08-15

set -euo pipefail

# Clean directory before compiling
rm -f *.o *.so


# To print the beta coefficients at the end of each cycle, add "-DWRITE_BETAS" to PKG_CFLAGS
# "-Wall -Wextra" flags compile with warnings enabled

# Compile ENCE and MICE DURR C code into a single .so
PKG_CFLAGS="-fopenmp -O2 -march=native -funroll-loops -flto -Wall -Wextra" \
PKG_LIBS="-fopenmp -flto" \
R CMD SHLIB elastic_net_functions.c elastic_net.c ENCE_impute.c MICE_DURR_impute.c ENCE_MICE_R_wrapper.c \
  -o ENCE_MICE_DURR_impute.so

# Clean up .o files
rm -f *.o

echo "Compiled ENCE_MICE_DURR_impute.so"