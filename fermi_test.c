#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "fermi.h"

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
    int n = 8;
    int G_order[3] = {0, 1, 2};
    int G_neg[3] = {1, 1, 1};

    double num_electrons = 1.0;
    double E_Fermi = 0.0;
    int err = FindFermi(n, num_bands, num_electrons, G_order, G_neg, Efn, &E_Fermi);
    printf("Got E_Fermi = %f\n", E_Fermi);
    if (err != CTETRA_BISECT_OK) {
        printf("Error = %d returned from FindFermi.\n", err);
        return err;
    }

    double expected = 13.0;
    if (E_Fermi < 12.0 || E_Fermi > 14.0) {
        printf("Incorrect E_Fermi; got %f, expected %f\n", E_Fermi, expected);
        return 1;
    }

    num_electrons = 0.5;
    err = FindFermi(n, num_bands, num_electrons, G_order, G_neg, Efn, &E_Fermi);
    printf("Got E_Fermi = %f\n", E_Fermi);
    if (err != CTETRA_BISECT_OK) {
        printf("Error = %d returned from FindFermi.\n", err);
        return err;
    }

    double tol = 1e-8;
    expected = 6.0;
    if (fabs(E_Fermi - expected) > tol) {
        printf("Incorrect E_Fermi; got %f, expected %f\n", E_Fermi, expected);
        return 1;
    }

    printf("Fermi test passed.\n");
    return 0;
}
