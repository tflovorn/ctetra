#ifndef CTETRA_SUM_H
#define CTETRA_SUM_H

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "input.h"
#include "ecache.h"
#include "submesh.h"
#include "fermi.h"
#include "weights.h"

double SumEnergy(double *E_Fermi, InputFn Efn, int n, int num_bands, double num_electrons, gsl_matrix *R, bool use_cache);

#endif // CTETRA_SUM_H
