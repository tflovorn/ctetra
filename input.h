#ifndef CTETRA_INPUT_H
#define CTETRA_INPUT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef void (*InputFn)(double k[3], gsl_vector *values);

typedef void (*UEInputFn)(double k[3], gsl_vector *evalues, gsl_matrix_complex *evecs);

#endif // CTETRA_INPUT_H
