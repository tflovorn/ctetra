#include "sum.h"

double SumEnergy(double *E_Fermi, InputFn Efn, int n, int num_bands, double num_electrons, gsl_matrix *R, bool use_cache) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);

    EnergyCache *Ecache = init_EnergyCache(n, num_bands, G_order, G_neg, Efn, use_cache);

    int err = FindFermi(n, num_bands, num_electrons, G_order, G_neg, Ecache, E_Fermi);
    if (err != CTETRA_BISECT_OK) {
        printf("Error: FindFermi failed with error code = %d\n", err);
        free_EnergyCache(Ecache)
        exit(EXIT_FAILURE);
    }
    free_EnergyCache(Ecache)
    return 0.0;
}
