#include "submesh.h"

void submesh_ijk_to_k(int n, int i, int j, int k, double k_opt[3]) {
    double step = 1.0 / ((double)n);
    k_opt[0] = ((double)i) * step;
    k_opt[1] = ((double)j) * step;
    k_opt[2] = ((double)k) * step;
}

void get_k_orig(double k_opt[3], int G_order[3], int G_neg[3], double k_orig[3]) {
    int i;
    for (i = 0; i < 3; i++) {
        k_orig[i] = k_opt[G_order[i]] * ((double)G_neg[i]);
    }
}
