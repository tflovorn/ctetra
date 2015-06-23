#ifndef CTETRA_SUBMESH_H
#define CTETRA_SUBMESH_H

void submesh_ijk_to_k(int n, int i, int j, int k, double k_opt[3]);

void get_k_orig(double k_opt[3], int G_order[3], int G_neg[3], double k_orig[3]);

#endif // CTETRA_SUBMESH_H
