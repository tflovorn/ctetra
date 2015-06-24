#include "weights.h"

void WeightsAtK(double E_Fermi, int i, int j, int k, EnergyCache *Ecache, double *result) {
    int n = Ecache->n;
    double num_tetra = (double)(6*n*n*n);
    int num_bands = Ecache->num_bands;
    double* c = (double*)malloc(num_bands * sizeof(double));
    int band_index;
    for (band_index = 0; band_index < num_bands; band_index++) {
        c[band_index] = 0.0;
    }
    double contrib, y, t;

    int tetras[6][4] = {{1, 2, 3, 6}, {1, 3, 5, 6}, {3, 5, 6, 7}, {3, 6, 7, 8}, {3, 4, 6, 8}, {2, 3, 4, 6}};
    int subcell_points[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
    int subcell_num = 0;
    int **subcells = subcells_around_ijk(n, i, j, k, &subcell_num);

    int subcell_index, tetra_index, vertex_index, point_index;
    int sci, scj, sck, pi, pj, pk;
    bool ijk_seen;

    gsl_vector **tetra_Es = (gsl_vector**)malloc(4 * sizeof(gsl_vector*));
    int **tetra_ks = (int**)malloc(4 * sizeof(int*));
    double *sorted_Es = (double*)malloc(4 * sizeof(double));
    int **sorted_ks = (int**)malloc(4 * sizeof(int*));
    double *ws = (double*)malloc(4 * sizeof(double));
    for (vertex_index = 0; vertex_index < 4; vertex_index++) {
        tetra_Es[vertex_index] = gsl_vector_alloc(num_bands);
        tetra_ks[vertex_index] = (int*)malloc(3 * sizeof(int));
        sorted_ks[vertex_index] = (int*)malloc(3 * sizeof(int));
    }

    // Iterate over subcells neighboring (i, j, k).
    for (subcell_index = 0; subcell_index < subcell_num; subcell_index++) {
        sci = subcells[subcell_index][0];
        scj = subcells[subcell_index][1];
        sck = subcells[subcell_index][2];
        // Iterate over the tetrahedra covering this subcell.
        for (tetra_index = 0; tetra_index < 6; tetra_index++) {
            // Skip tetrahedra that do not include (i, j, k) as a vertex.
            ijk_seen = false;
            for (vertex_index = 0; vertex_index < 4; vertex_index++) {
                point_index = tetras[tetra_index][vertex_index] - 1;
                pi = sci + subcell_points[point_index][0];
                pj = scj + subcell_points[point_index][1];
                pk = sck + subcell_points[point_index][2];
                if (pi == i && pj == j && pk == k) {
                    ijk_seen = true;
                }
            }
            if (!ijk_seen) {
                continue;
            }
            // If we get here, want to include this tetrahedron.
            // Collect band energies at the vertices.
            for (vertex_index = 0; vertex_index < 4; vertex_index++) {
                point_index = tetras[tetra_index][vertex_index] - 1;
                pi = sci + subcell_points[point_index][0];
                pj = scj + subcell_points[point_index][1];
                pk = sck + subcell_points[point_index][2];
                energy_from_cache(Ecache, pi, pj, pk, tetra_Es[vertex_index]);
                tetra_ks[vertex_index][0] = pi;
                tetra_ks[vertex_index][1] = pj;
                tetra_ks[vertex_index][2] = pk;
            }
            // For each band, sort vertices by energy and calculate the
            // associated weight.
            for (band_index = 0; band_index < num_bands; band_index++) {
                // Initialize sorting arrays for this band.
                for (vertex_index = 0; vertex_index < 4; vertex_index++) {
                    sorted_Es[vertex_index] = gsl_vector_get(tetra_Es[vertex_index], band_index);
                    sorted_ks[vertex_index][0] = tetra_ks[vertex_index][0];
                    sorted_ks[vertex_index][1] = tetra_ks[vertex_index][1];
                    sorted_ks[vertex_index][2] = tetra_ks[vertex_index][2];
                }
                // Sort Es and sort ks by associated Es.
                sortEsKs(sorted_Es, sorted_ks);
                // Get the weight contributions for this tetrahedron's vertices.
                WeightContrib(E_Fermi, sorted_Es[0], sorted_Es[1], sorted_Es[2], sorted_Es[3], num_tetra, ws);
                // Add the weight contribution from vertex (i, j, k) to the total.
                for (vertex_index = 0; vertex_index < 4; vertex_index++) {
                    pi = sorted_ks[vertex_index][0];
                    pj = sorted_ks[vertex_index][1];
                    pk = sorted_ks[vertex_index][2];
                    if (pi == i && pj == j && pk == k) {
                        // Shouldn't get here more than once per band_index.
                        // TODO - add check to verify?
                        contrib = ws[vertex_index];
                        y = contrib - c[band_index];
                        t = result[band_index] + y;
                        c[band_index] = (t - result[band_index]) - y;
                        result[band_index] = t;
                    }
                }
            }
        }
    }

    for (vertex_index = 0; vertex_index < 4; vertex_index++) {
        free(sorted_ks[vertex_index]);
        free(tetra_ks[vertex_index]);
        gsl_vector_free(tetra_Es[vertex_index]);
    }
    free(ws);
    free(sorted_ks);
    free(sorted_Es);
    free(tetra_ks);
    free(tetra_Es);

    for (subcell_index = 0; subcell_index < subcell_num; subcell_index++) {
        free(subcells[subcell_index]);
    }
    free(subcells);
    free(c);
}

// Assumes that Es and ks are length-4 arrays.
void sortEsKs(double *Es, int **ks) {
    // insertion sort
    int i, j;
    double tmpE;
    int *tmpk;
	for (i = 1; i < 4; i++) {
		for (j = i; j > 0 && Es[j-1] > Es[j]; j--) {
			// swap j, j-1
            tmpE = Es[j-1];
			Es[j-1] = Es[j];
            Es[j] = tmpE;

            tmpk = ks[j-1];
            ks[j-1] = ks[j];
            ks[j] = tmpk;
		}
	}
}

// Return the specified tetrahedron's contribution to the integration
// weights at the k-points of the tetrahedron's vertices; i.e. return
// a list with elements w_{bandIndex, kN, tetra}, where the elements of the
// returned list range over kN values in the order specified by tetra.
// The calculation of the tetrahedron contribution to w_{nj} is implemented
// as described in BJA94 Appendix B and Section V.
//
// E_Fermi = Fermi energy of the system.
//
// TODO docs
void WeightContrib(double E_Fermi, double E1, double E2, double E3, double E4, double num_tetra, double *ws) {
	double qnt = 1.0 / (4.0 * num_tetra);
	if (E_Fermi <= E1) {
		int i;
		for (i = 0; i < 4; i++) {
			ws[i] = 0.0;
		}
	} else if (E1 <= E_Fermi && E_Fermi <= E2) {
        if (E1 == E2) {
            int i;
            for (i = 0; i < 4; i++) {
                ws[i] = 0.0;
            }
        }
		double C = qnt * (E_Fermi - E1) * (E_Fermi - E1) * (E_Fermi - E1) / ((E2 - E1) * (E3 - E1) * (E4 - E1));
		ws[0] = C * (4.0 - (E_Fermi-E1)*(1.0/(E2-E1)+1.0/(E3-E1)+1.0/(E4-E1)));
		ws[1] = C * (E_Fermi - E1) / (E2 - E1);
		ws[2] = C * (E_Fermi - E1) / (E3 - E1);
		ws[3] = C * (E_Fermi - E1) / (E4 - E1);
	} else if (E2 <= E_Fermi && E_Fermi <= E3) {
        if (E2 == E3) {
            double C1 = qnt * (E_Fermi - E1) * (E_Fermi - E1) / ((E4 - E1) * (E3 - E1));
            ws[0] = C1 + C1*(E4 - E_Fermi)/(E4 - E1);
            ws[1] = C1;
            ws[2] = C1*(E_Fermi - E1)/(E3 - E1);
            ws[3] = C1*(E_Fermi - E1)/(E4 - E1);
        }
		double C1 = qnt * (E_Fermi - E1) * (E_Fermi - E1) / ((E4 - E1) * (E3 - E1));

		double C2_num = qnt * (E_Fermi - E1) * (E_Fermi - E2) * (E3 - E_Fermi);
		double C2_denom = (E4 - E1) * (E3 - E2) * (E3 - E1);
		double C2 = C2_num / C2_denom;

		double C3_num = qnt * (E_Fermi - E2) * (E_Fermi - E2) * (E4 - E_Fermi);
		double C3_denom = (E4 - E2) * (E3 - E2) * (E4 - E1);
		double C3 = C3_num / C3_denom;

		ws[0] = C1 + (C1+C2)*(E3-E_Fermi)/(E3-E1) + (C1+C2+C3)*(E4-E_Fermi)/(E4-E1);
		ws[1] = C1 + C2 + C3 + (C2+C3)*(E3-E_Fermi)/(E3-E2) + C3*(E4-E_Fermi)/(E4-E2);
		ws[2] = (C1+C2)*(E_Fermi-E1)/(E3-E1) + (C2+C3)*(E_Fermi-E2)/(E3-E2);
		ws[3] = (C1+C2+C3)*(E_Fermi-E1)/(E4-E1) + C3*(E_Fermi-E2)/(E4-E2);
	} else if (E3 <= E_Fermi && E_Fermi <= E4) {
        if (E3 == E4) {
            int i;
            for (i = 0; i < 4; i++) {
                ws[i] = qnt;
            }
        }
		double C = qnt * (E4 - E_Fermi) * (E4 - E_Fermi) * (E4 - E_Fermi) / ((E4 - E1) * (E4 - E2) * (E4 - E3));
		ws[0] = qnt - C*(E4-E_Fermi)/(E4-E1);
		ws[1] = qnt - C*(E4-E_Fermi)/(E4-E2);
		ws[2] = qnt - C*(E4-E_Fermi)/(E4-E3);
		ws[3] = qnt - C*(4.0-(1.0/(E4-E1)+1.0/(E4-E2)+1.0/(E4-E3))*(E4-E_Fermi));
	} else {
		// E_Fermi >= E4
		int i;
		for (i = 0; i < 4; i++) {
			ws[i] = qnt;
		}
	}

	addCurvatureCorrection(E_Fermi, E1, E2, E3, E4, num_tetra, ws);
}

// Return a list of the curvature corrections to the k-point weight
// contributions from the specified tetrahedron. The band energies at the
// vertices of the tetrahedron are given in sorted order by Es; the returned
// corrections are given in the same order.
// TODO fix doc
void addCurvatureCorrection(double E_Fermi, double E1, double E2, double E3, double E4, double num_tetra, double *ws) {
	double D_T = DosContrib(E_Fermi, E1, E2, E3, E4, num_tetra);
	double sumEs = E1 + E2 + E3 + E4;
	ws[0] += (D_T / 40.0) * (sumEs - 4.0*E1);
	ws[1] += (D_T / 40.0) * (sumEs - 4.0*E2);
	ws[2] += (D_T / 40.0) * (sumEs - 4.0*E3);
	ws[3] += (D_T / 40.0) * (sumEs - 4.0*E4);
}
