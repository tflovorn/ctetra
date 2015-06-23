#include "submesh.h"

void submesh_ijk_to_k(int n, int i, int j, int k, double k_opt[3]) {
    double step = 1.0 / ((double)n);
    k_opt[0] = ((double)i) * step;
    k_opt[1] = ((double)j) * step;
    k_opt[2] = ((double)k) * step;
}

void get_k_orig(double k_opt[3], int G_order[3], int G_neg[3], double k_orig[3]) {
    int i;
    for (i = 0; i < 3; i++) {
        k_orig[i] = k_opt[G_order[i]] * ((double)G_neg[i]);
    }
}

void MinMaxVals(int n, int num_bands, int G_order[3], int G_neg[3], InputFn Efn, double *emin, double *emax) {
    double minval = 0.0;
    double maxval = 0.0;
    double this_min, this_max;
    double k_opt[3] = {0.0, 0.0, 0.0};
    double k_orig[3] = {0.0, 0.0, 0.0};
    gsl_vector *energies = gsl_vector_alloc(num_bands);
    int i, j, k;
    for (k = 0; k < n+1; k++) {
        for (j = 0; j < n+1; j++) {
            for (i = 0; i < n+1; i++) {
                submesh_ijk_to_k(n, i, j, k, k_opt);
                get_k_orig(k_opt, G_order, G_neg, k_orig);
                Efn(k_orig, energies);
                gsl_vector_minmax(energies, &this_min, &this_max);
                if (i == 0 && j == 0 && k == 0) {
                    minval = this_min;
                    maxval = this_max;
                } else {
                    if (this_min < minval) {
                        minval = this_min;
                    }
                    if (this_max > maxval) {
                        maxval = this_max;
                    }
                }
            }
        }
    }
    *emin = minval;
    *emax = maxval;
}
