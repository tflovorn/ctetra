#ifndef CTETRA_FERMI_H
#define CTETRA_FERMI_H

#include <time.h>
#include <math.h>
#include "submesh.h"
#include "numstates.h"

#define CTETRA_BISECT_OK 0
#define CTETRA_BISECT_BRACKET 1
#define CTETRA_BISECT_MAXITER 2

typedef double (*BisectFn)(double x);

int FindFermi(int n, int num_bands, double num_electrons, int G_order[3], int G_neg[3], EnergyCache *Ecache, double *E_Fermi);

int Bisect(double *result, BisectFn fn, double a, double b, double toly, int maxiter);

#endif // CTETRA_FERMI_H
