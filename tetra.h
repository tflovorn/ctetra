#ifndef CTETRA_TETRA_H
#define CTETRA_TETRA_H

#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "ecache.h"

typedef double (*tetra_SumFn)(double E, double E1, double E2, double E3, double E4, double num_tetra);

double tetra_SumTetra(tetra_SumFn F, double E, EnergyCache *Ecache);

void sortEs(double Es[4]);

#endif // CTETRA_TETRA_H
