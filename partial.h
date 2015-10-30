#ifndef CTETRA_PARTIAL_H
#define CTETRA_PARTIAL_H

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include "input.h"
#include "evcache.h"
#include "submesh.h"

double Gauss_PartialDos_Eig(double E, double sigma, int orig_index, int eig_index, EvecCache *evCache);

double Gauss_PartialDos(double E, double sigma, int orig_index, EvecCache *evCache);

double** Gauss_PartialDosList(UEInputFn UEfn, int n, double sigma, int num_bands, gsl_matrix *R, double **Es, int num_dos);

#endif //CTETRA_PARTIAL_H
