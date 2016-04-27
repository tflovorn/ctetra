#ifndef CTETRA_NUMSTATES_H
#define CTETRA_NUMSTATES_H

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "input.h"
#include "submesh.h"
#include "ecache.h"
#include "tetra.h"

double NumStates(double E, EnergyCache *Ecache);

double NumStatesContrib(double E, double E1, double E2, double E3, double E4, double num_tetra);

#endif // CTETRA_NUMSTATES_H
