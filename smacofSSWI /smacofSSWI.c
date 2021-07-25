#include "../include/smacof.h"

int main() {
  double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double w[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double x[8] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  int n = 4, p = 2, itmax = 100;
  bool verbose = true, speedup = false;
  double eps = 1e-15;
  (void)smacofSSWI(delta, w, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  double *dist = (double *)calloc((size_t)6, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  (void)smacofPrintMatrix(1, 6, 15, 10, delta);
  (void)smacofPrintMatrix(1, 6, 15, 10, dist);
  (void)smacofGradientW(dist, w, delta, n, p, x);
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  free(dist);
  return EXIT_SUCCESS;
}

void smacofSSWI(double *delta, const double *w, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, smid = 0.0, snew = 0.0;
  // normalize delta
  (void)smacofNormDeltaW(w, delta, m);
  // Compute the MP inverse of V and test for irreducibility
  double *vinv = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofInvertVW(w, vinv, n, m);
  // compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  // scale initial configuration
  (void)smacofScaleXW(x, dist, w, delta, m, np);
  // compute initial stress
  sold = smacofLossW(dist, w, delta, m);
  while (true > false) {
    (void)smacofGuttmanW(dist, w, delta, vinv, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    smid = smacofLossW(dist, w, delta, m);
    snew = smacofLossW(dist, w, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f smid = %15.10f snew = %15.10f\n", itel,
             sold, smid, snew);
    }
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(dist);
  free(vinv);
  return;
}


void smacofIntervalW(const double *dist, const double *w, double *delta,
                     int m) {
  double acum = 0.0, bcum = 0.0, ccum = 0.0, dcum = 0.0, wcum = 0.0;
  double s = 0.0, deltamin = INFINITY;
  for (int k = 0; k < m; k++) {
    acum += w[k] * delta[k];
    bcum += w[k] * dist[k];
    wcum += w[k];
    deltamin = MIN(deltamin, delta[k]);
  }
  acum /= wcum;
  bcum /= wcum;
  for (int k = 0; k < m; k++) {
    ccum += w[k] * (delta[k] - acum) * (dist[k] - bcum);
    dcum += w[k] * SQUARE(delta[k] - acum);
  }
  double a0 = ccum / dcum;
  double b0 = -(a0 * acum - bcum);
  double c0 = a0 * deltamin + b0;
  if ((a0 > 0.0) && ((a0 * deltamin + b0) > 0.0)) {
    printf("case 1: a = %15.10f b = %15.10f min = %15.10f\n", a0, b0, c0);
    s = 0.0;
    for (int k = 0; k < m; k++) {
      delta[k] = a0 * delta[k] + b0;
      s += w[k] * SQUARE(delta[k]);
    }
    double r = sqrt(((double)m) / s);
    for (int k = 0; k < m; k++) {
      delta[k] *= r;
    }
    return;
  } else {
    ccum = 0.0;
    dcum = 0.0;
    for (int k = 0; k < m; k++) {
      ccum += w[k] * (delta[k] - deltamin) * dist[k];
      dcum += w[k] * SQUARE(delta[k] - deltamin);
    }
    double a1 = ccum / dcum;
    double b1 = -a1 * deltamin;
    double c1 = a1 * deltamin + b1;
    printf("case 2: a = %15.10f b = %15.10f min = %15.10f\n", a1, b1, c1);
    s = 0.0;
    for (int k = 0; k < m; k++) {
      delta[k] = a1 * delta[k] + b1;
      s += w[k] * SQUARE(delta[k]);
    }
    double r = sqrt(((double)m) / s);
    for (int k = 0; k < m; k++) {
      delta[k] *= r;
    }
    return;
  }
}
