#ifndef CWANNIER_EVCACHE_H
#define CWANNIER_EVCACHE_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "input.h"
#include "ecache.h"
#include "submesh.h"

typedef struct {
    int na;
    int nb;
    int nc;
    int num_bands;
    int G_order[3];
    int G_neg[3];
    UEInputFn UEfn;
    gsl_vector **energies;
    gsl_matrix_complex **evecs;
} EvecCache;

EvecCache* init_EvecCache(int na, int nb, int nc, int num_bands, int G_order[3], int G_neg[3], UEInputFn UEfn);

void free_EvecCache(EvecCache *evCache);

void EvecCache_MinMaxVals(EvecCache *evCache, double *emin, double *emax);

EnergyCache* copy_to_Ecache(EvecCache *evCache);

#endif //CWANNIER_EVCACHE_H
