#include "numstates.h"

// Return n(E), the total number of states with energy <= E summed over
// all tetrahedra and band indices.
// The calculation of n(E) is implemented as described in BJA94 Appendix A.
//
// TODO doc
// Uses Kahan summation for improved accuracy on dense mesh.
double NumStates(double E, int n, int num_bands, int G_order[3], int G_neg[3], InputFn Efn) {
    double num_tetra = (double)(6*n*n*n);
    double k_opt[3] = {0.0, 0.0, 0.0};
    double k_orig[3] = {0.0, 0.0, 0.0};
    double result = 0.0;
    double c = 0.0;
    double contrib, y, t;
    int i, j, k, point_index, tetra_index, pi, pj, pk, band_index, vertex_index;
    int subcell_points[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
    int tetras[6][4] = {{1, 2, 3, 6}, {1, 3, 5, 6}, {3, 5, 6, 7}, {3, 6, 7, 8}, {3, 4, 6, 8}, {2, 3, 4, 6}};
    double tetraEs[4] = {0.0, 0.0, 0.0, 0.0};

    gsl_vector **energies = malloc(8 * sizeof(gsl_vector*));
    for (i = 0; i < 8; i++) {
        energies[i] = gsl_vector_alloc(num_bands);
        gsl_vector_set_zero(energies[i]);
    }

    // Iterate over tetrahedra.
    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                // Collect band energies at vertices.
                for (point_index = 0; point_index < 8; point_index++) {
                    pi = i + subcell_points[point_index][0];
                    pj = j + subcell_points[point_index][1];
                    pk = k + subcell_points[point_index][2];
                    submesh_ijk_to_k(n, pi, pj, pk, k_opt);
                    get_k_orig(k_opt, G_order, G_neg, k_orig);
                    // TODO replace this with cached val if required.
                    Efn(k_orig, energies[point_index]);
                }
                // Calculate contributions for each tetrahedron.
                for (tetra_index = 0; tetra_index < 6; tetra_index++) {
                    // For each band: sort values, then add contributions.
                    for (band_index = 0; band_index < num_bands; band_index++) {
                        for (vertex_index = 0; vertex_index < 4; vertex_index++) {
                            point_index = tetras[tetra_index][vertex_index] - 1;
                            tetraEs[vertex_index] = gsl_vector_get(energies[point_index], band_index);
                        }
                        sortEs(tetraEs);
                        contrib = NumStatesContrib(E, tetraEs[0], tetraEs[1], tetraEs[2], tetraEs[3], num_tetra);
                        y = contrib - c;
                        t = result + y;
                        c = (t - result) - y;
                        result = t;
                    }
                }
            }
        }
    }

    for (i = 0; i < 8; i++) {
        gsl_vector_free(energies[i]);
    }
    free(energies);

    return result;
}

void sortEs(double Es[4]) {
    // insertion sort
    int i, j;
    double tmp;
	for (i = 1; i < 4; i++) {
		for (j = i; j > 0 && Es[j-1] > Es[j]; j--) {
			// swap j, j-1
            tmp = Es[j-1];
			Es[j-1] = Es[j];
            Es[j] = tmp;
		}
	}
}

// Return the contribution to the number of states with energy less than
// or equal to E (i.e. the integrated density of states n(E)) from the
// (tetrahedron, band index) pair with the given energies at the vertices.
// The calculation of n(E) is implemented as described in BJA94 Appendix A.
//
// E1, E2, E3, E4 = energies at the vertices of the tetrahedron, in ascending
// order.
//
// num_tetra = total number of tetrahedra in the full Brillouin zone.
// Equal to (volume of tetrahedron) / (volume of full BZ).
double NumStatesContrib(double E, double E1, double E2, double E3, double E4, double num_tetra) {
	if (E <= E1) {
		return 0.0;
	} else if (E1 < E && E < E2) {
		return (1.0 / num_tetra) * (E - E1) * (E - E1) * (E - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1));
	} else if (E2 < E && E < E3) {
		double fac = (1.0 / num_tetra) / ((E3 - E1) * (E4 - E1));
		double esq = (E2-E1)*(E2-E1) + 3.0*(E2-E1)*(E-E2) + 3.0*(E-E2)*(E-E2);
		double ecub = -(((E3 - E1) + (E4 - E2)) / ((E3 - E2) * (E4 - E2))) * (E - E2) * (E - E2) * (E - E2);
		return fac * (esq + ecub);
	} else if (E3 < E && E < E4) {
		return (1.0 / num_tetra) * (1.0 - (E4-E)*(E4-E)*(E4-E)/((E4-E1)*(E4-E2)*(E4-E3)));
	} else {
		// E >= E4
		return (1.0 / num_tetra);
	}
}
