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

void OptimizeGs(gsl_matrix *R, int G_order[3], int G_neg[3]) {
    int permutations[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    int signs[8][3] = {{1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1}, {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};
    double opt_36 = 0.0;
    double k3_to_k6 = 0.0;
    int pi, si, i;
    gsl_matrix *this_R = gsl_matrix_alloc(3, 3);
    gsl_vector *Ri = gsl_vector_alloc(3);
    gsl_vector *mR0 = gsl_vector_alloc(3);
    gsl_vector *R1 = gsl_vector_alloc(3);
    gsl_vector *mR2 = gsl_vector_alloc(3);
    for (pi = 0; pi < 6; pi++) {
        for (si = 0; si < 8; si++) {
            // Get next (permutation, sign) pair.
            for (i = 0; i < 3; i++) {
                // this_R[perm[i], :] = sign[i] * R[i, :]
                // TODO check gsl error values?
                gsl_matrix_get_row(Ri, R, i);
                gsl_vector_scale(Ri, signs[si][i]);
                gsl_matrix_set_row(this_R, permutations[pi][i], Ri);
            }
            // Calculate 3->6 distance.
            // mR0 = -1 * this_R[0, :]
            gsl_matrix_get_row(mR0, this_R, 0);
            gsl_vector_scale(mR0, -1.0);
            // R1 = this_R[1, :]
            gsl_matrix_get_row(R1, this_R, 1);
            // mR2 = -1 * this_R[2, :]
            gsl_matrix_get_row(mR2, this_R, 2);
            gsl_vector_scale(mR2, -1.0);
            // k3_to_k6 = norm(mR0 + R1 + mR2)
            gsl_vector_add(mR0, R1);
            gsl_vector_add(mR0, mR2);
            k3_to_k6 = gsl_blas_dnrm2(mR0);
            //printf("got k3_to_k6 = %f\n", k3_to_k6);
            // Is this 3->6 distance smaller?
            if ((pi == 0 && si == 0) || k3_to_k6 < opt_36) {
                opt_36 = k3_to_k6;
                for (i = 0; i < 3; i++) {
                    G_order[i] = permutations[pi][i];
                    G_neg[i] = signs[si][i];
                }
            }
        }
    }
    gsl_matrix_free(this_R);
    gsl_vector_free(Ri);
    gsl_vector_free(mR0);
    gsl_vector_free(R1);
    gsl_vector_free(mR2);
    //printf("finished optimizeGs; k3_to_k6 = %f\n", k3_to_k6);
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
