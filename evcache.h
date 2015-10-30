#ifndef CWANNIER_EVCACHE_H
#define CWANNIER_EVCACHE_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "input.h"
#include "ecache.h"
#include "submesh.h"

typedef struct {
    int n;
    int num_bands;
    int G_order[3];
    int G_neg[3];
    UEInputFn UEfn;
    gsl_vector **energies;
    gsl_matrix_complex **evecs;
} EvecCache;

EvecCache* init_EvecCache(int n, int num_bands, int G_order[3], int G_neg[3], UEInputFn UEfn);

void free_EvecCache(EvecCache *evCache);

void EvecCache_MinMaxVals(int n, int num_bands, EvecCache *evCache, double *emin, double *emax);

#endif //CWANNIER_EVCACHE_H
