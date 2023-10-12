#include "../include/smacof.h"

int main() {
  double delta[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double x[8] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
  int n = 4, p = 2, itmax = 1000;
  bool verbose = true, speedup = false;
  double eps = 1e-15;
  (void)smacofSSUI(delta, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  double *dist = (double *)calloc((size_t)6, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  (void)smacofPrintMatrix(1, 6, 15, 10, delta);
  (void)smacofPrintMatrix(1, 6, 15, 10, dist);
  (void)smacofGradientU(dist, delta, n, p, x);
  (void)smacofPrintMatrix(4, 2, 15, 10, x);
  free(dist);
  return EXIT_SUCCESS;
}

void smacofSSUI(double *delta, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, smid = 0.0, snew = 0.0;
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
    smid = smacofLossU(dist, delta, m);
    (void)smacofIntervalU(dist, delta, m);
    snew = smacofLossU(dist, delta, m);
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
  return;
}


void smacofIntervalU(const double *dist, double *delta, int m) {
  double acum = 0.0, bcum = 0.0, ccum = 0.0, dcum = 0.0;
  double s = 0.0, deltamin = INFINITY;
  for (int k = 0; k < m; k++) {
    acum += delta[k];
    bcum += dist[k];
    deltamin = MIN(deltamin, delta[k]);
  }
//  acum /= (double)m;
//  bcum /= (double)m;
  for (int k = 0; k < m; k++) {
    ccum += (delta[k] - acum) * (dist[k] - bcum);
    dcum += SQUARE(delta[k] - acum);
  }
  double a0 = ccum / dcum;
  double b0 = -(a0 * acum - bcum);
  if ((a0 > 1e-15) && ((a0 * deltamin + b0) > 1e-15)) {
    s = 0.0;
    for (int k = 0; k < m; k++) {
      delta[k] = a0 * delta[k] + b0;
      s += SQUARE(delta[k]);
    }
  } else {
    ccum = 0.0;
    dcum = 0.0;
    for (int k = 0; k < m; k++) {
      ccum += (delta[k] - deltamin) * dist[k];
      dcum += SQUARE(delta[k] - deltamin);
    }
    double a1 = ccum / dcum;
    double b1 = -a1 * deltamin;
    if (a1 > 1e-15) {
      s = 0.0;
      for (int k = 0; k < m; k++) {
        delta[k] = a0 * delta[k] + b0;
        s += SQUARE(delta[k]);
      }
    } else {
      return;
    }
  }
  double r = sqrt(1 / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}
