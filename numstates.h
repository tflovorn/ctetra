#ifndef CTETRA_NUMSTATES_H
#define CTETRA_NUMSTATES_H

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "submesh.h"
#include "sum.h"

double NumStates(double E, int n, int num_bands, int G_order[3], int G_neg[3], InputFn Efn);

void sortEs(double Es[4]);

double NumStatesContrib(double E, double E1, double E2, double E3, double E4, double num_tetra);

#endif // CTETRA_NUMSTATES_H
