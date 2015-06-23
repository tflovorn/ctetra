#ifndef CTETRA_SUBMESH_H
#define CTETRA_SUBMESH_H

#include "sum.h"

void submesh_ijk_to_k(int n, int i, int j, int k, double k_opt[3]);

void get_k_orig(double k_opt[3], int G_order[3], int G_neg[3], double k_orig[3]);

void MinMaxVals(int n, int num_bands, int G_order[3], int G_neg[3], InputFn Efn, double *emin, double *emax);

#endif // CTETRA_SUBMESH_H
