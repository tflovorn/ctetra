#include "weights.h"

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
