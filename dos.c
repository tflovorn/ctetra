#include "dos.h"

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
