#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include "numstates.h"

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
    int G_order[3] = {0, 1, 2};
    int G_neg[3] = {1, 1, 1};
    bool use_cache = true;

    EnergyCache *Ecache = init_EnergyCache(na, nb, nc, num_bands, G_order, G_neg, Efn, use_cache);

    double E = 0.0;
    double count = NumStates(E, Ecache);
    double eps = 1e-9;
    double expected_E0 = 0.0;
    if (fabs(count - expected_E0) > eps) {
        printf("Incorrect occupation; got %f, expected %f\n", count, expected_E0);
        return 1;
    }

    E = 6.0;
    count = NumStates(E, Ecache);
    double expected_E6 = 0.5;
    if (fabs(count - expected_E6) > eps) {
        printf("Incorrect occupation; got %f, expected %f\n", count, expected_E6);
        return 1;
    }

    E = 12.0;
    count = NumStates(E, Ecache);
    double expected_E12 = 1.0;
    if (fabs(count - expected_E12) > eps) {
        printf("Incorrect occupation; got %f, expected %f\n", count, expected_E12);
        return 1;
    }

    E = 26.0;
    count = NumStates(E, Ecache);
    double expected_E24 = 2.0;
    if (fabs(count - expected_E24) > eps) {
        printf("Incorrect occupation; got %f, expected %f\n", count, expected_E24);
        return 1;
    }

    printf("Numstates test passed.\n");
    return 0;
}
