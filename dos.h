#ifndef CTETRA_DOS_H
#define CTETRA_DOS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "input.h"
#include "ecache.h"
#include "tetra.h"
#include "submesh.h"
#include "fermi.h"

double* Tetra_AllDosList(InputFn Efn, int na, int nb, int nc, int num_bands, gsl_matrix *R, double *Es, int num_dos);
double* Tetra_DosList(InputFn Efn, int na, int nb, int nc, int num_bands, gsl_matrix *R, double *Es, int num_dos);
double* Tetra_DosEnergyDerivList(InputFn Efn, int na, int nb, int nc, int num_bands, gsl_matrix *R, double *Es, int num_dos, double num_electrons, double *fermi, double *dos_fermi, double *dos_deriv_fermi);
double Tetra_TotalDos(double E, EnergyCache *Ecache);
double Tetra_TotalDosEnergyDeriv(double E, EnergyCache *Ecache);
double DosContrib(double E, double E1, double E2, double E3, double E4, double num_tetra);
double DosEnergyDerivContrib(double E, double E1, double E2, double E3, double E4, double num_tetra);

#endif // CTETRA_DOS_H
