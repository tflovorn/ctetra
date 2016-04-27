#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "sum.h"

int main(int argc, char *argv[]) {
    int num_bands = 2;
    double t = 1.0;
    double E0 = 6.0;
    double deltaE = 14.0;
    void Efn(double k[3], gsl_vector *energies) {
        double cx = cos(2.0 * M_PI * k[0]);
        double cy = cos(2.0 * M_PI * k[1]);
        double cz = cos(2.0 * M_PI * k[2]);
        double tk = -2.0 * t * (cx + cy + cz);
        int i;
        double E0_band;
        for (i = 0; i < num_bands; i++) {
            E0_band = E0 + ((double)i) * deltaE;
            gsl_vector_set(energies, i, E0_band + tk);
        }
    }
    int na = 8;
    int nb = 8;
    int nc = 8;
    bool use_cache = true;

    gsl_matrix *R = gsl_matrix_calloc(3, 3); // R -> all zeros
    // Overall scale for R doesn't matter.
    gsl_matrix_set(R, 0, 0, 1.0);
    gsl_matrix_set(R, 1, 1, 1.0);
    gsl_matrix_set(R, 2, 2, 1.0);

    double num_electrons = 1.0;
    double E_Fermi = 0.0;
    double energy = SumEnergy(&E_Fermi, Efn, na, nb, nc, num_bands, num_electrons, R, use_cache);
    printf("Got E_Fermi = %f\n", E_Fermi);
    printf("Got energy = %f\n", energy);

    double tol = 1e-8;
    double expected = 6.0;
    if (fabs(energy - expected) > tol) {
        printf("Incorrect energy; got %f, expected %f\n", energy, expected);
        return 1;
    }

    printf("Sum test passed.\n");
    return 0;
}
