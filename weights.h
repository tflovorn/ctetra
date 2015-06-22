#ifndef CTETRA_WEIGHTS_H
#define CTETRA_WEIGHTS_H

#include "dos.h"

void WeightContrib(double E_Fermi, double E1, double E2, double E3, double E4, double num_tetra, double *ws);

void addCurvatureCorrection(double E_Fermi, double E1, double E2, double E3, double E4, double num_tetra, double *ws);

#endif // CTETRA_WEIGHTS_H
