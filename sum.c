#include "sum.h"

double SumEnergy(double *E_Fermi, InputFn Efn, int n, int num_bands, double num_electrons, gsl_matrix *R, bool use_cache) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);
    int err = FindFermi(n, num_bands, num_electrons, G_order, G_neg, Efn, E_Fermi);
    if (err != CTETRA_BISECT_OK) {
        printf("Error: FindFermi failed with error code = %d\n", err);
        exit(EXIT_FAILURE);
    }
    return 0.0;
}
