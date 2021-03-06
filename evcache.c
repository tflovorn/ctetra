#include "evcache.h"

EvecCache* init_EvecCache(int na, int nb, int nc, int num_bands, int G_order[3], int G_neg[3], UEInputFn UEfn) {
    int num_ks = (na+1)*(nb+1)*(nc+1);
    int i, j, k;
    EvecCache *evCache = (EvecCache*)malloc(sizeof(EvecCache));
    evCache->na = na;
    evCache->nb = nb;
    evCache->nc = nc;
    evCache->num_bands = num_bands;
    for (i = 0; i < 3; i++) {
        evCache->G_order[i] = G_order[i];
        evCache->G_neg[i] = G_neg[i];
    }
    evCache->UEfn = UEfn;

    // Initialize and calculate the value of the band energies and 
    // eigenvectors at each k-point.
    gsl_vector **energies = (gsl_vector**)malloc(num_ks * sizeof(gsl_vector*));
    evCache->energies = energies;
    gsl_matrix_complex **evecs = (gsl_matrix_complex**)malloc(num_ks * sizeof(gsl_matrix_complex*));
    evCache->evecs = evecs;

    for (k = 0; k < nc+1; k++) {
        for (j = 0; j < nb+1; j++) {
            for (i = 0; i < na+1; i++) {
                int kN = submesh_ijk_index(na, nb, nc, i, j, k);
                evCache->energies[kN] = gsl_vector_calloc(num_bands);
                evCache->evecs[kN] = gsl_matrix_complex_calloc(num_bands, num_bands);

                double k_opt[3] = {0, 0, 0};
                double k_orig[3] = {0, 0, 0};
                submesh_ijk_to_k(na, nb, nc, i, j, k, k_opt);
                get_k_orig(k_opt, G_order, G_neg, k_orig);
                UEfn(k_orig, evCache->energies[kN], evCache->evecs[kN]);
                verify_sort(evCache->energies[kN]);
            }
        }
    }

    return evCache;
}

void free_EvecCache(EvecCache *evCache) {
    int na = evCache->na;
    int nb = evCache->nb;
    int nc = evCache->nc;
    int num_ks = (na+1)*(nb+1)*(nc+1);
    int kN;
    for (kN = 0; kN < num_ks; kN++) {
        gsl_vector_free(evCache->energies[kN]);
        gsl_matrix_complex_free(evCache->evecs[kN]);
    }
    free(evCache->energies);
    free(evCache->evecs);
    free(evCache);
}

void EvecCache_MinMaxVals(EvecCache *evCache, double *emin, double *emax) {
    int na = evCache->na;
    int nb = evCache->nb;
    int nc = evCache->nc;
    double minval = 0.0;
    double maxval = 0.0;
    double this_min, this_max;
    gsl_vector *energies;

    int i, j, k;
    for (k = 0; k < nc+1; k++) {
        for (j = 0; j < nb+1; j++) {
            for (i = 0; i < na+1; i++) {
                int kN = submesh_ijk_index(na, nb, nc, i, j, k);
                energies = evCache->energies[kN];

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

EnergyCache* copy_to_Ecache(EvecCache *evCache) {
    EnergyCache *Ecache = malloc(sizeof(EnergyCache));
    Ecache->na = evCache->na;
    Ecache->nb = evCache->nb;
    Ecache->nc = evCache->nc;
    Ecache->num_bands = evCache->num_bands;
    Ecache->energies = evCache->energies;

    Ecache->use_cache = true;
    Ecache->Efn = NULL;

    int i;
    for (i = 0; i < 3; i++) {
        Ecache->G_order[i] = 0;
        Ecache->G_neg[i] = 0;
    }

    return Ecache;
}
