#include "../include/smacof.h"

void smacofNormW(double *w, const int m) {
  double s = 0.0; 
  for (int k = 0; k < m; k++) {
    s += w[k];
  }
  for (int k = 0; k < m; k++) {
    w[k] /= s;
  }
  return;
}

void smacofNormDeltaW(const double *w, double *delta, const int m) {
  double s = 0.0; 
  for (int k = 0; k < m; k++) {
    s += w[k] * SQUARE(delta[k]);
  }
  double r = sqrt(1.0 / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}

void smacofInvertVW(const double *w, double *vinv, const int n, const int m) {
  double ni = 1.0 / (double)n, s = 0.0;
  for (int k = 0; k < m; k++) {
    vinv[k] = -w[k] + ni;
  }
  double *vsum = (double *)calloc((size_t)n, sizeof(double));
  for (int i = 1; i <= n; i++) {
    s = 0.0;
    for (int j = i + 1; j <= n; j++) {
      s += w[SINDEX(j, i, n)];
    }
    for (int j = 1; j < i; j++) {
      s += w[SINDEX(i, j, n)];
    }
    vsum[VINDEX(i)] = s + ni;
  }
  double *vkrw = (double *)calloc((size_t)n, sizeof(double));
  for (int k = 1; k <= n; k++) {
    for (int j = 1; j <= n; j++) {
      if (j == k) {
        vkrw[VINDEX(k)] = 0.0;
      }
      if (j < k) {
        vkrw[VINDEX(j)] = vinv[SINDEX(k, j, n)];
      }
      if (j > k) {
        vkrw[VINDEX(j)] = vinv[SINDEX(j, k, n)];
      }
    }
    s = vsum[VINDEX(k)];
    if (s < 1e-15) {
      printf("Error: W is not irreducible.\n");
      exit(EXIT_FAILURE);
    }
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= i; j++) {
        if (i == j) {
          vsum[VINDEX(i)] -= SQUARE(vkrw[VINDEX(i)]) / s;
        } else {
          vinv[SINDEX(i, j, n)] -= vkrw[VINDEX(i)] * vkrw[VINDEX(j)] / s;
        }
      }
    }
    for (int j = 1; j < k; j++) {
      vinv[SINDEX(k, j, n)] /= s;
    }
    for (int i = k + 1; i <= n; i++) {
      vinv[SINDEX(i, k, n)] /= s;
    }
    vsum[VINDEX(k)] = -1.0 / s;
  }
  for (int k = 0; k < m; k++) {
    vinv[k] = -vinv[k] - ni;
  }
  free(vsum);
  free(vkrw);
  return;
}

void smacofScaleXW(double *x, double *dist, const double *w,
                   const double *delta, const int m, const int np) {
  double sd1 = 0.0, sd2 = 0.0;
  for (int k = 0; k < m; k++) {
    sd1 += w[k] * dist[k] * delta[k];
    sd2 += w[k] * SQUARE(dist[k]);
  }
  double lbd = sd1 / sd2;
  for (int k = 0; k < m; k++) {
    dist[k] *= lbd;
  }
  for (int k = 0; k < np; k++) {
    x[k] *= lbd;
  }
}

double smacofLossW(const double *dist, const double *w, const double *delta,
                   const int m) {
  double stress = 0.0;
  for (int k = 0; k < m; k++) {
    stress += w[k] * SQUARE(delta[k] - dist[k]);
  }
  return sqrt(stress);
}

void smacofGuttmanW(const double *dist, const double *w, const double *delta,
                    const double *vinv, const int n, const int p,
                    const bool speedup, double *x) {
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
        y[VINDEX(i)] += w[k] * delta[k] *
                        (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]) / dist[k];
      }
    }
    double *z = calloc((size_t)n, sizeof(double));
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
        z[VINDEX(i)] += vinv[k] * (y[VINDEX(j)] - y[VINDEX(i)]);
      }
      if (speedup) {
        z[VINDEX(i)] = 2 * z[VINDEX(i)] - x[MINDEX(i, s, n)];
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = z[VINDEX(i)];
    }
    free(y);
    free(z);
  }
}

void smacofGradientW(double *dist, double *w, double *delta, const int n,
                       const int p, double *x) { 
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
        y[VINDEX(i)] += w[k] * (1 - delta[k] / dist[k]) *
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

