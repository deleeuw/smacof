#include "../include/smacof.h"

int main() {
  double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double x[8] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  int n = 4, p = 2, itmax = 100;
  bool verbose = true, speedup = false;
  double eps = 1e-15;
  (void)smacofSSUR(delta, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  double *dist = (double *)calloc((size_t)6, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  (void)smacofGradientU(dist, delta, n, p, x);
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  free(dist);
  return EXIT_SUCCESS;
}

void smacofSSUR(double *delta, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, snew = 0.0;
  // normalize delta
  (void)smacofNormDeltaU(delta, m);
  // compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  // scale initial configuration
  (void)smacofScaleXU(x, dist, delta, m, np);
  // compute initial stress
  sold = smacofLossU(dist, delta, m);
  while (true > false) {
    (void)smacofGuttmanU(dist, delta, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    snew = smacofLossU(dist, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f snew = %15.10f\n", itel, sold, snew);
    }
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(dist);
  return;
}

