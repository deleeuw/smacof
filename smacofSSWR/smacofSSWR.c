#include "../include/smacof.h"

int main() {
  double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double w[6] = {1.0, 1.0, 0.0, 2.0, 2.0, 2.0};
  double x[8] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  int n = 4, p = 2, itmax = 100;
  bool verbose = true, speedup = false;
  double eps = 1e-10;
  (void)smacofSSWR(delta, w, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  double *dist = (double *)calloc((size_t) 6, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  (void)smacofGradientW(dist, w, delta, n, p, x); 
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  free(dist);
  return EXIT_SUCCESS;
}

// Single Matrix, Symmetric, Weighted Stress, Ratio Transform

void smacofSSWR(double *delta, double *w, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, snew = 0.0;
  // normalize w
  (void)smacofNormW(w, m);
  // normalize delta
  (void)smacofNormDeltaW(w, delta, m);
  // Compute the MP inverse of V and test for irreducibility
  double *vinv = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofInvertVW(w, vinv, n, m);
  // Compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  // Scale initial configuration and distances
  (void)smacofScaleXW(x, dist, w, delta, m, np);
  // Compute initial stress
  sold = smacofLossW(dist, w, delta, m);
  while (true) {
    (void)smacofGuttmanW(dist, w, delta, vinv, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    snew = smacofLossW(dist, w, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f snew = %15.10f\n", itel, sold, snew);
    }
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(vinv);
  free(dist);
  return;
}

