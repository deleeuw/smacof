#include "../include/smacof.h"

void smacofNormDeltaU(double *delta, const int m) {
  double s = 0.0;
  for (int k = 0; k < m; k++) {
    s += SQUARE(delta[k]);
  }
  double r = sqrt(1.0 / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}

void smacofScaleXU(double *x, double *dist, const double *delta, const int m,
                   const int np) {
  double sd1 = 0.0, sd2 = 0.0;
  for (int k = 0; k < m; k++) {
    sd1 += dist[k] * delta[k];
    sd2 += SQUARE(dist[k]);
  }
  double lbd = sd1 / sd2;
  for (int k = 0; k < m; k++) {
    dist[k] *= lbd;
  }
  for (int k = 0; k < np; k++) {
    x[k] *= lbd;
  }
}

double smacofLossU(const double *dist, const double *delta, const int m) {
  double stress = 0.0;
  for (int k = 0; k < m; k++) {
    stress += SQUARE(delta[k] - dist[k]);
  }
  return stress;
}

void smacofGuttmanU(const double *dist, const double *delta, const int n,
                    const int p, const bool speedup, double *x) {
  int k;
  for (int s = 1; s <= p; s++) {
    double *y = calloc((size_t)n, sizeof(double));
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        if (j == i) {
          continue;
        }
        if (i > j) {
          k = SINDEX(i, j, n);
        }
        if (j > i) {
          k = SINDEX(j, i, n);
        }
        if (dist[k] < 1e-15) {
          continue;
        }
        y[VINDEX(i)] +=
            delta[k] * (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]) / dist[k];
      }
      y[VINDEX(i)] /= (double)n;
      if (speedup) {
        y[VINDEX(i)] = 2 * y[VINDEX(i)] - x[MINDEX(i, s, n)];
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = y[VINDEX(i)];
    }
    free(y);
  }
  return;
}

void smacofGradientU(double *dist, double *delta, const int n, const int p,
                     double *x) {
  for (int s = 1; s <= p; s++) {
    double *y = (double *)calloc((size_t)n, sizeof(double));
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        int k = 0;
        if (j == i) {
          continue;
        }
        if (i > j) {
          k = SINDEX(i, j, n);
        }
        if (j > i) {
          k = SINDEX(j, i, n);
        }
        if (dist[k] < 1e-15) {
          continue;
        }
        y[VINDEX(i)] += (1 - delta[k] / dist[k]) *
                        (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = y[VINDEX(i)];
    }
    free(y);
  }
  return;
}
