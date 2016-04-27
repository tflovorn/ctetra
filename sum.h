#ifndef CTETRA_SUM_H
#define CTETRA_SUM_H

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "input.h"
#include "ecache.h"
#include "evcache.h"
#include "submesh.h"
#include "fermi.h"
#include "weights.h"

double SumEnergy(double *E_Fermi, InputFn Efn, int na, int nb, int nc, int num_bands, double num_electrons, gsl_matrix *R, bool use_cache);

double SumEnergyFixedFermi(double E_Fermi, InputFn Efn, int na, int nb, int nc, int num_bands, gsl_matrix *R, bool use_cache);

double with_cache_SumEnergyFixedFermi(double E_Fermi, EnergyCache *Ecache);

double** partial_num_states(UEInputFn UEfn, int na, int nb, int nc, int num_bands, double num_total_electrons, gsl_matrix *R, double **Es, int num_E, double *E_Fermi, double **num_states_Fermi);

double one_band_n(double E_Fermi, int band_index, EvecCache *evCache);

#endif // CTETRA_SUM_H
