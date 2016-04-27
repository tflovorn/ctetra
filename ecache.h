#ifndef CWANNIER_ECACHE_H
#define CWANNIER_ECACHE_H

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "input.h"

typedef struct {
    int na;
    int nb;
    int nc;
    int num_bands;
    int G_order[3];
    int G_neg[3];
    InputFn Efn;
    bool use_cache;
    gsl_vector **energies;
} EnergyCache;

#include "submesh.h"

EnergyCache* init_EnergyCache(int na, int nb, int nc, int num_bands, int G_order[3], int G_neg[3], InputFn Efn, bool use_cache);

void free_EnergyCache(EnergyCache *Ecache);

void energy_from_cache(EnergyCache *Ecache, int i, int j, int k, gsl_vector *energies);

void verify_sort(gsl_vector *energies);

#endif // CWANNIER_ECACHE_H
