#ifndef CTETRA_SUBMESH_H
#define CTETRA_SUBMESH_H

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "ecache.h"

void submesh_ijk_to_k(int n, int i, int j, int k, double k_opt[3]);

int submesh_ijk_index(int n, int i, int j, int k);

void get_k_orig(double k_opt[3], int G_order[3], int G_neg[3], double k_orig[3]);

void MinMaxVals(int n, int num_bands, int G_order[3], int G_neg[3], EnergyCache *Ecache, double *emin, double *emax);

void OptimizeGs(gsl_matrix *R, int G_order[3], int G_neg[3]);

#endif // CTETRA_SUBMESH_H
