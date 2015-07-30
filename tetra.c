#include "tetra.h"

// Return \sum_{k} F(E), the value of F summed over all tetrahedra.
//
// TODO doc
// Uses Kahan summation for improved accuracy on dense mesh.
double tetra_SumTetra(tetra_SumFn F, double E, int n, int num_bands, EnergyCache *Ecache) {
    double num_tetra = (double)(6*n*n*n);
    double result = 0.0;
    double c = 0.0;
    double contrib, y, t;
    int i, j, k, point_index, tetra_index, pi, pj, pk, band_index, vertex_index;
    int subcell_points[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
    int tetras[6][4] = {{1, 2, 3, 6}, {1, 3, 5, 6}, {3, 5, 6, 7}, {3, 6, 7, 8}, {3, 4, 6, 8}, {2, 3, 4, 6}};
    double tetraEs[4] = {0.0, 0.0, 0.0, 0.0};

    gsl_vector **energies = malloc(8 * sizeof(gsl_vector*));
    for (i = 0; i < 8; i++) {
        energies[i] = gsl_vector_calloc(num_bands);
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
                    energy_from_cache(Ecache, pi, pj, pk, energies[point_index]);
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
                        contrib = F(E, tetraEs[0], tetraEs[1], tetraEs[2], tetraEs[3], num_tetra);
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


