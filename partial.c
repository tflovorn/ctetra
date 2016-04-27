#include "partial.h"

// Compute partial density of states corresponding to band index `orig_index`.
// Compute only the contribution from eigenvalue index `eig_index`.
// Perform calculation using Gaussian broadening.
double Gauss_PartialDos_Eig(double E, double sigma, int orig_index, int eig_index, EvecCache *evCache) {
    int na = evCache->na;
    int nb = evCache->nb;
    int nc = evCache->nc;
    int Nk = na*nb*nc; // not (n+1)^3 since we will avoid the repeated points
    double fac = sqrt(2.0*M_PI) * sigma * Nk;
    double Ek;
    gsl_matrix_complex *U;
    gsl_complex U_part;
    int i, j, k, k_index;
    double weight, gauss;
    double pdos = 0.0, this_pdos = 0.0;

    for (k = 0; k < nc; k++) {
        for (j = 0; j < nb; j++) {
            for (i = 0; i < na; i++) {
                k_index = submesh_ijk_index(na, nb, nc, i, j, k);

                U = evCache->evecs[k_index];
                U_part = gsl_matrix_complex_get(U, orig_index, eig_index);
                weight = gsl_complex_abs2(U_part);

                Ek = gsl_vector_get(evCache->energies[k_index], eig_index);
                gauss = exp(-pow(E - Ek, 2.0) / (2.0 * pow(sigma, 2.0)));

                this_pdos = weight * gauss;
                pdos += this_pdos;
            }
        }
    }

    return pdos / fac;
}

// Compute partial density of states corresponding to band index `orig_index`.
// Include contriubutions from all eigenvalue indices.
// Perform calculation using Gaussian broadening.
double Gauss_PartialDos(double E, double sigma, int orig_index, EvecCache *evCache) {
    int eig_index;
    int num_bands = evCache->num_bands;
    double pdos = 0.0, this_pdos = 0.0;

    for (eig_index = 0; eig_index < num_bands; eig_index++) {
        this_pdos = Gauss_PartialDos_Eig(E, sigma, orig_index, eig_index, evCache);
        pdos += this_pdos;
    }

    return pdos;
}

// Return a list of lists of partial density of states values between
// the minimum and maximum energy eigenvalues. The first index in the
// returned list corresponds to the state index in the original basis;
// the number of these entries is num_bands. The second index in the returned
// list has length num_dos, corresponding to equally spaced energy values
// between the minimum and maximum energy eigenvalues.
// The energy values used are stored in Es; after this function returns, *Es
// is a list with length num_dos.
// NOTE this behavior is different than Es in dos.c -- TODO make consistent?
double** Gauss_PartialDosList(UEInputFn UEfn, int na, int nb, int nc, double sigma, int num_bands, gsl_matrix *R, double **Es, int num_dos) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    int band_index, E_index;
    double emin, emax;

    OptimizeGs(R, G_order, G_neg);
    EvecCache *evCache = init_EvecCache(na, nb, nc, num_bands, G_order, G_neg, UEfn);

    EvecCache_MinMaxVals(evCache, &emin, &emax);
    double step = (emax - emin) / (num_dos - 1);

    double **dos_vals = malloc(num_bands * sizeof(double*));
    for (band_index = 0; band_index < num_bands; band_index++) {
        dos_vals[band_index] = malloc(num_dos * sizeof(double));
    }
    *Es = malloc(num_dos * sizeof(double));
    
    for (band_index = 0; band_index < num_bands; band_index++) {
        for (E_index = 0; E_index < num_dos; E_index++) {
            double E = emin + E_index*step;
            (*Es)[E_index] = E;
            dos_vals[band_index][E_index] = Gauss_PartialDos(E, sigma, band_index, evCache);
        }

    }
    return dos_vals;
}
