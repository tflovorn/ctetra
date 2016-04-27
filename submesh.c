#include "submesh.h"

void submesh_ijk_to_k(int na, int nb, int nc, int i, int j, int k, double k_opt[3]) {
    double step_a = 1.0 / ((double)na);
    double step_b = 1.0 / ((double)nb);
    double step_c = 1.0 / ((double)nc);
    k_opt[0] = ((double)i) * step_a;
    k_opt[1] = ((double)j) * step_b;
    k_opt[2] = ((double)k) * step_c;
}

int submesh_ijk_index(int na, int nb, int nc, int i, int j, int k) {
    return i + j*(na+1) + k*(na+1)*(nb+1);
}

int **subcells_around_ijk(int na, int nb, int nc, int i, int j, int k, int *subcell_num) {
    *subcell_num = 0;
	// For each subcell, point (i, j, k) has the number in BJA94 Fig. 5
	// corresponding to the subcell number.
    // First pass: get number of subcells that will be returned.
	// Subcell #1
	if (i != na && j != nb && k != nc) {
        (*subcell_num)++;
	}
	// Subcell #2
	if (i != 0 && j != nb && k != nc) {
        (*subcell_num)++;
	}
	// Subcell #3
	if (i != na && j != 0 && k != nc) {
        (*subcell_num)++;
	}
	// Subcell #4
	if (i != 0 && j != 0 && k != nc) {
        (*subcell_num)++;
	}
	// Subcell #5
	if (i != na && j != nb && k != 0) {
        (*subcell_num)++;
	}
	// Subcell #6
	if (i != 0 && j != nb && k != 0) {
        (*subcell_num)++;
	}
	// Subcell #7
	if (i != na && j != 0 && k != 0) {
        (*subcell_num)++;
	}
	// Subcell #8
	if (i != 0 && j != 0 && k != 0) {
        (*subcell_num)++;
	}
    // Allocate return value.
    int **subcells = (int**)malloc((*subcell_num) * sizeof(int*));
    int subcell_index = 0;
    // Collect subcells to return.
	// Subcell #1
	if (i != na && j != nb && k != nc) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i;
        subcells[subcell_index][1] = j;
        subcells[subcell_index][2] = k;
        subcell_index++;
	}
	// Subcell #2
	if (i != 0 && j != nb && k != nc) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i - 1;
        subcells[subcell_index][1] = j;
        subcells[subcell_index][2] = k;
        subcell_index++;
	}
	// Subcell #3
	if (i != na && j != 0 && k != nc) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i;
        subcells[subcell_index][1] = j - 1;
        subcells[subcell_index][2] = k;
        subcell_index++;
	}
	// Subcell #4
	if (i != 0 && j != 0 && k != nc) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i - 1;
        subcells[subcell_index][1] = j - 1;
        subcells[subcell_index][2] = k;
        subcell_index++;
	}
	// Subcell #5
	if (i != na && j != nb && k != 0) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i;
        subcells[subcell_index][1] = j;
        subcells[subcell_index][2] = k - 1;
        subcell_index++;
	}
	// Subcell #6
	if (i != 0 && j != nb && k != 0) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i - 1;
        subcells[subcell_index][1] = j;
        subcells[subcell_index][2] = k - 1;
        subcell_index++;
	}
	// Subcell #7
	if (i != na && j != 0 && k != 0) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i;
        subcells[subcell_index][1] = j - 1;
        subcells[subcell_index][2] = k - 1;
        subcell_index++;
	}
	// Subcell #8
	if (i != 0 && j != 0 && k != 0) {
        subcells[subcell_index] = (int*)malloc(3*sizeof(int));
        subcells[subcell_index][0] = i - 1;
        subcells[subcell_index][1] = j - 1;
        subcells[subcell_index][2] = k - 1;
        subcell_index++;
	}
    return subcells;
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

void MinMaxVals(EnergyCache *Ecache, double *emin, double *emax) {
    int na = Ecache->na;
    int nb = Ecache->nb;
    int nc = Ecache->nc;
    int num_bands = Ecache->num_bands;
    double minval = 0.0;
    double maxval = 0.0;
    double this_min, this_max;
    gsl_vector *energies = gsl_vector_alloc(num_bands);
    int i, j, k;
    for (k = 0; k < nc+1; k++) {
        for (j = 0; j < nb+1; j++) {
            for (i = 0; i < na+1; i++) {
                energy_from_cache(Ecache, i, j, k, energies);
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
    gsl_vector_free(energies);
}
