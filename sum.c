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
