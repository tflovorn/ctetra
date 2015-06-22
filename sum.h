#ifndef CTETRA_SUM_H
#define CTETRA_SUM_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>

typedef void (*InputFn)(double k[3], gsl_vector *values);

double SumEnergy(double *E_Fermi, InputFn Efn, int n, double num_electrons, double R[3][3], bool use_cache);

#endif // CTETRA_SUM_H
