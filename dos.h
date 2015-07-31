#ifndef CTETRA_DOS_H
#define CTETRA_DOS_H

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "input.h"
#include "ecache.h"
#include "tetra.h"
#include "submesh.h"

double* Tetra_DosList(InputFn Efn, int n, int num_bands, gsl_matrix *R, double *Es, int num_dos);
double Tetra_TotalDos(double E, EnergyCache *Ecache, int n, int num_bands);
double DosContrib(double E, double E1, double E2, double E3, double E4, double num_tetra);

#endif // CTETRA_DOS_H
