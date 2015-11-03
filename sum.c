#include "sum.h"

double SumEnergy(double *E_Fermi, InputFn Efn, int n, int num_bands, double num_electrons, gsl_matrix *R, bool use_cache) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);

    EnergyCache *Ecache = init_EnergyCache(n, num_bands, G_order, G_neg, Efn, use_cache);

    int err = FindFermi(n, num_bands, num_electrons, Ecache, E_Fermi);
    if (err != CTETRA_BISECT_OK) {
        printf("Error: FindFermi failed with error code = %d\n", err);
        free_EnergyCache(Ecache);
        exit(EXIT_FAILURE);
    }

    double energy = with_cache_SumEnergyFixedFermi(*E_Fermi, Ecache, n, num_bands);

    free_EnergyCache(Ecache);

    return energy;
}

double SumEnergyFixedFermi(double E_Fermi, InputFn Efn, int n, int num_bands, gsl_matrix *R, bool use_cache) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);

    EnergyCache *Ecache = init_EnergyCache(n, num_bands, G_order, G_neg, Efn, use_cache);

    double energy = with_cache_SumEnergyFixedFermi(E_Fermi, Ecache, n, num_bands);

    free_EnergyCache(Ecache);

    return energy;
}

double with_cache_SumEnergyFixedFermi(double E_Fermi, EnergyCache *Ecache, int n, int num_bands) {
    double result = 0.0;
    double c = 0.0;
    double contrib, y, t;
    int i, j, k, band_index;
    double *this_ws = (double*)malloc(num_bands * sizeof(double));
    gsl_vector *this_Es = gsl_vector_calloc(num_bands);
    // Iterate over all k-points and collect their contributions to
    // the energy.
    // The equivalent indices 0 and n are both included since the
    // weights for points with these indices are generated from distinct
    // places in the Brillouin zone.
    for (k = 0; k < n+1; k++) {
        for (j = 0; j < n+1; j++) {
            for (i = 0; i < n+1; i++) {
                for (band_index = 0; band_index < num_bands; band_index++) {
                    this_ws[band_index] = 0.0;
                }
                WeightsAtK(E_Fermi, i, j, k, Ecache, this_ws);
                energy_from_cache(Ecache, i, j, k, this_Es);
                for (band_index = 0; band_index < num_bands; band_index++) {
                    contrib = this_ws[band_index] * gsl_vector_get(this_Es, band_index);
                    y = contrib - c;
                    t = result + y;
                    c = (t - result) - y;
                    result = t;
                }
            }
        }
    }
    gsl_vector_free(this_Es);
    free(this_ws);

    return result;
}

double** partial_num_states(UEInputFn UEfn, int n, int num_bands, double num_total_electrons, gsl_matrix *R, double **Es, int num_E, double *E_Fermi, double **num_states_Fermi) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);

    EvecCache *evCache = init_EvecCache(n, num_bands, G_order, G_neg, UEfn);
    EnergyCache *Ecache = copy_to_Ecache(evCache);

    int err = FindFermi(n, num_bands, num_total_electrons, Ecache, E_Fermi);
    if (err != CTETRA_BISECT_OK) {
        printf("Error: FindFermi failed with error code = %d\n", err);
        free_EnergyCache(Ecache);
        exit(EXIT_FAILURE);
    }

    double emin, emax;
    EvecCache_MinMaxVals(n, num_bands, evCache, &emin, &emax);
    double step = (emax - emin) / (num_E - 1);

    double **num_states_vals = malloc(num_bands * sizeof(double*));
    int band_index, E_index;
    for (band_index = 0; band_index < num_bands; band_index++) {
        num_states_vals[band_index] = malloc(num_E * sizeof(double));
    }
    *Es = malloc(num_E * sizeof(double));

    *num_states_Fermi = malloc(num_bands * sizeof(double));
    for (band_index = 0; band_index < num_bands; band_index++) {
        (*num_states_Fermi)[band_index] = one_band_n(*E_Fermi, band_index, evCache);
    }

    for (band_index = 0; band_index < num_bands; band_index++) {
        for (E_index = 0; E_index < num_E; E_index++) {
            double E = emin + E_index*step;
            (*Es)[E_index] = E;
            num_states_vals[band_index][E_index] = one_band_n(E, band_index, evCache);
        }
    }

    free(Ecache);
    free_EvecCache(evCache);

    return num_states_vals;
}

double one_band_n(double E_Fermi, int band_index, EvecCache *evCache) {
    double result = 0.0;
    double c = 0.0;
    double contrib, y, t;
    int i, j, k, eig_band_index, kN;
    int num_bands = evCache->num_bands;
    int n = evCache->n;
    double prob_val;
    double *this_ws = (double*)malloc(num_bands * sizeof(double));
    EnergyCache *Ecache = copy_to_Ecache(evCache);
    gsl_matrix_complex *U;
    // Iterate over all k-points and collect their contributions to
    // the number of states sum.
    // The equivalent indices 0 and n are both included since the
    // weights for points with these indices are generated from distinct
    // places in the Brillouin zone.
    for (k = 0; k < n+1; k++) {
        for (j = 0; j < n+1; j++) {
            for (i = 0; i < n+1; i++) {
                for (eig_band_index = 0; eig_band_index < num_bands; eig_band_index++) {
                    this_ws[eig_band_index] = 0.0;
                }
                WeightsAtK(E_Fermi, i, j, k, Ecache, this_ws);
                kN = submesh_ijk_index(evCache->n, i, j, k);
                U = evCache->evecs[kN];
                for (eig_band_index = 0; eig_band_index < num_bands; eig_band_index++) {
                    prob_val = gsl_complex_abs2(gsl_matrix_complex_get(U, band_index, eig_band_index));
                    contrib = this_ws[eig_band_index] * prob_val;
                    y = contrib - c;
                    t = result + y;
                    c = (t - result) - y;
                    result = t;
                }
            }
        }
    }
    free(Ecache);
    free(this_ws);

    return result;
}
