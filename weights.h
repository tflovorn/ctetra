#ifndef CTETRA_WEIGHTS_H
#define CTETRA_WEIGHTS_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "ecache.h"
#include "dos.h"

void WeightsAtK(double E_Fermi, int i, int j, int k, EnergyCache *Ecache, double *result);

void sortEsKs(double *Es, int **ks);

void WeightContrib(double E_Fermi, double E1, double E2, double E3, double E4, double num_tetra, double *ws);

void addCurvatureCorrection(double E_Fermi, double E1, double E2, double E3, double E4, double num_tetra, double *ws);

#endif // CTETRA_WEIGHTS_H
