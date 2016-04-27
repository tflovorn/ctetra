#include "dos.h"

// Return a list of density of states values between the minimum and maximum
// energy eigenvalues, which has length num_dos.
// The energy values used are stored in Es, which should have length num_dos.
double* Tetra_AllDosList(InputFn Efn, int na, int nb, int nc, int num_bands, gsl_matrix *R, double *Es, int num_dos) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);

    bool use_cache = true;
    EnergyCache *Ecache = init_EnergyCache(na, nb, nc, num_bands, G_order, G_neg, Efn, use_cache);

    double emin, emax;
    MinMaxVals(Ecache, &emin, &emax);
    double step = (emax - emin) / (num_dos - 1);

    double *dos_vals = malloc(num_dos * sizeof(double));

    int i = 0;
    for (i = 0; i < num_dos; i++) {
        double E = emin + i*step;
        Es[i] = E;
        dos_vals[i] = Tetra_TotalDos(E, Ecache);
    }
    return dos_vals;
}

// Return a list of density of states values at the energies given in Es,
// which has length num_dos.
double* Tetra_DosList(InputFn Efn, int na, int nb, int nc, int num_bands, gsl_matrix *R, double *Es, int num_dos) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);

    bool use_cache = true;
    EnergyCache *Ecache = init_EnergyCache(na, nb, nc, num_bands, G_order, G_neg, Efn, use_cache);

    double *dos_vals = malloc(num_dos * sizeof(double));

    int i = 0;
    for (i = 0; i < num_dos; i++) {
        double E = Es[i];
        dos_vals[i] = Tetra_TotalDos(E, Ecache);
    }
    return dos_vals;
}

// Return a list of density of states derivative values between the minimum and maximum
// energy eigenvalues, which has length num_dos.
// The energy values used are stored in Es, which should have length num_dos.
double* Tetra_DosEnergyDerivList(InputFn Efn, int na, int nb, int nc, int num_bands, gsl_matrix *R, double *Es, int num_dos, double num_electrons, double *fermi, double *dos_fermi, double *dos_deriv_fermi) {
    int G_order[3] = {0, 0, 0};
    int G_neg[3] = {0, 0, 0};
    OptimizeGs(R, G_order, G_neg);

    bool use_cache = true;
    EnergyCache *Ecache = init_EnergyCache(na, nb, nc, num_bands, G_order, G_neg, Efn, use_cache);

    int err = FindFermi(num_electrons, Ecache, fermi);
    if (err != CTETRA_BISECT_OK) {
        printf("Error: FindFermi failed with error code = %d\n", err);
        free_EnergyCache(Ecache);
        exit(EXIT_FAILURE);
    }
    *dos_fermi = Tetra_TotalDos(*fermi, Ecache);
    *dos_deriv_fermi = Tetra_TotalDosEnergyDeriv(*fermi, Ecache);

    double emin, emax;
    MinMaxVals(Ecache, &emin, &emax);
    double step = (emax - emin) / (num_dos - 1);

    double *dos_deriv_vals = malloc(num_dos * sizeof(double));

    int i = 0;
    for (i = 0; i < num_dos; i++) {
        double E = emin + i*step;
        Es[i] = E;
        dos_deriv_vals[i] = Tetra_TotalDosEnergyDeriv(E, Ecache);
    }
    return dos_deriv_vals;
}

// Return the density of states at energy E.
double Tetra_TotalDos(double E, EnergyCache *Ecache) {
    return tetra_SumTetra(DosContrib, E, Ecache);
}

// Return the energy derivative of the density of states at energy E.
double Tetra_TotalDosEnergyDeriv(double E, EnergyCache *Ecache) {
    return tetra_SumTetra(DosEnergyDerivContrib, E, Ecache);
}

// Return the contribution to the density of states at energy E from the
// (tetrahedron, band index) pair with the given energies at the vertices.
// The calculation of D_T(E) is implemented as described in BJA94 Appendix C.
//
// E1, E2, E3, E4 = energies at the vertices of the tetrahedron, in ascending
// order.
//
// num_tetra = total number of tetrahedra in the full Brillouin zone.
// Equal to (volume of tetrahedron) / (volume of full BZ).
double DosContrib(double E, double E1, double E2, double E3, double E4, double num_tetra) {
	if (E <= E1) {
		return 0.0;
	} else if (E1 <= E && E <= E2) {
        if (E1 == E2) {
            return 0.0;
        }
		return (1.0 / num_tetra) * 3.0 * (E - E1) * (E - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1));
	} else if (E2 <= E && E <= E3) {
        if (E2 == E3) {
            return (1.0 / num_tetra) * 3.0 * (E2 - E1) / ((E3 - E1) * (E4 - E1));
        }
		double fac = (1.0 / num_tetra) / ((E3 - E1) * (E4 - E1));
		double elin = 3.0*(E2-E1) + 6.0*(E-E2);
		double esq = -3.0 * (((E3 - E1) + (E4 - E2)) / ((E3 - E2) * (E4 - E2))) * (E - E2) * (E - E2);
		return fac * (elin + esq);
	} else if (E3 <= E && E <= E4) {
        if (E3 == E4) {
            return 0.0;
        }
		return (1.0 / num_tetra) * 3.0 * (E4 - E) * (E4 - E) / ((E4 - E1) * (E4 - E2) * (E4 - E3));
	} else {
		// E >= E4
		return 0.0;
	}
}

// Return the contribution to the derivative of the density of states at with
// respect to energy at energy E from the (tetrahedron, band index) pair with
// the given energies at the vertices.
// The calculation of D_T(E) is implemented as described in BJA94 Appendix C.
//
// E1, E2, E3, E4 = energies at the vertices of the tetrahedron, in ascending
// order.
//
// num_tetra = total number of tetrahedra in the full Brillouin zone.
// Equal to (volume of tetrahedron) / (volume of full BZ).
double DosEnergyDerivContrib(double E, double E1, double E2, double E3, double E4, double num_tetra) {
	if (E <= E1) {
		return 0.0;
	} else if (E1 <= E && E <= E2) {
        if (E1 == E2) {
            return (1.0 / num_tetra) * 6.0 / ((E3 - E1) * (E4 - E1));
        }
		return (1.0 / num_tetra) * 6.0 * (E - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1));
	} else if (E2 <= E && E <= E3) {
		double fac = (1.0 / num_tetra) / ((E3 - E1) * (E4 - E1));
        if (E2 == E3) {
            return fac * (6.0 - 6.0*((E3 - E1) + (E4 - E2))/(E4 - E2));
        }
        return fac * (6.0 - 6.0*((E3 - E1) + (E4 - E2))*(E - E2)/((E3 - E2)*(E4 - E2)));
	} else if (E3 <= E && E <= E4) {
        if (E3 == E4) {
            return -(1.0 / num_tetra) * 6.0 / ((E4 - E1) * (E4 - E2));
        }
        return -(1.0 / num_tetra) * 6.0 * (E4 - E) / ((E4 - E1) * (E4 - E2) * (E4 - E3));
	} else {
		// E >= E4
		return 0.0;
	}
}
