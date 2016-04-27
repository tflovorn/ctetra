#include "ecache.h"

EnergyCache* init_EnergyCache(int na, int nb, int nc, int num_bands, int G_order[3], int G_neg[3], InputFn Efn, bool use_cache) {
    int num_ks = (na+1)*(nb+1)*(nc+1);
    int i, j, k;
    EnergyCache *Ecache = (EnergyCache*)malloc(sizeof(EnergyCache));
    Ecache->na = na;
    Ecache->nb = nb;
    Ecache->nc = nc;
    Ecache->num_bands = num_bands;
    for (i = 0; i < 3; i++) {
        Ecache->G_order[i] = G_order[i];
        Ecache->G_neg[i] = G_neg[i];
    }
    Ecache->Efn = Efn;
    Ecache->use_cache = use_cache;

    if (use_cache) {
        // Initialize and calculate the value of the band energies at
        // each k-point.
        gsl_vector **energies = (gsl_vector**)malloc(num_ks * sizeof(gsl_vector*));
        Ecache->energies = energies;

        for (k = 0; k < nc+1; k++) {
            for (j = 0; j < nb+1; j++) {
                for (i = 0; i < na+1; i++) {
                    int kN = submesh_ijk_index(Ecache->na, Ecache->nb, Ecache->nc, i, j, k);
                    Ecache->energies[kN] = gsl_vector_calloc(num_bands);

                    double k_opt[3] = {0, 0, 0};
                    double k_orig[3] = {0, 0, 0};
                    submesh_ijk_to_k(na, nb, nc, i, j, k, k_opt);
                    get_k_orig(k_opt, G_order, G_neg, k_orig);
                    Efn(k_orig, Ecache->energies[kN]);
                    verify_sort(Ecache->energies[kN]);
                }
            }
        }
    }
    return Ecache;
}

void free_EnergyCache(EnergyCache *Ecache) {
    if (Ecache->use_cache) {
        int na = Ecache->na;
        int nb = Ecache->nb;
        int nc = Ecache->nc;
        int num_ks = (na+1)*(nb+1)*(nc+1);
        int kN;
        for (kN = 0; kN < num_ks; kN++) {
            gsl_vector_free(Ecache->energies[kN]);
        }
        free(Ecache->energies);
    }
    free(Ecache);
}

void energy_from_cache(EnergyCache *Ecache, int i, int j, int k, gsl_vector *energies) {
    if (Ecache->use_cache) {
        int kN = submesh_ijk_index(Ecache->na, Ecache->nb, Ecache->nc, i, j, k);
        gsl_vector_memcpy(energies, Ecache->energies[kN]);
    } else {
        double k_opt[3] = {0, 0, 0};
        double k_orig[3] = {0, 0, 0};
        submesh_ijk_to_k(Ecache->na, Ecache->nb, Ecache->nc, i, j, k, k_opt);
        get_k_orig(k_opt, Ecache->G_order, Ecache->G_neg, k_orig);
        Ecache->Efn(k_orig, energies);
        verify_sort(energies);
    }
}

void verify_sort(gsl_vector *energies) {
    int i;
    int len = energies->size;
    for (i = 1; i < len; i++) {
        double at_i = gsl_vector_get(energies, i);
        double at_i_m1 = gsl_vector_get(energies, i-1);
        if (at_i_m1 > at_i) {
            printf("Error: energy eigenvalue vector not sorted.\n");
            exit(EXIT_FAILURE);
        }
    }
}
