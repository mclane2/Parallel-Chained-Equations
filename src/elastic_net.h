#ifndef ELASTIC_NET_H
#define ELASTIC_NET_H

int elastic_net_cov(const double *Z_in, const double *z_in, int n, int p, const double *lambdas, int n_lambda, double alpha,
                          double thresh, int maxit, double *beta_out, double *intercept_out, double *C, unsigned char *has_row);
#endif