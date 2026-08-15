#ifndef ELASTIC_NET_FUNCTIONS_H
#define ELASTIC_NET_FUNCTIONS_H

int standardise(double *Z, double *z, int n, int p, double *Zm, double *Zs, double *zm, double *zs);

double back_transform(double *beta, const double *Zm, const double *Zs, double zm, double zs, int p);

int coord_descent_cov(const double *Z, double *beta, int n, int p, double lambda, double alpha, double threshold, int maxit, int *is_active, int *active_idx,
                      double *g, double *C, char *has_row);

#endif