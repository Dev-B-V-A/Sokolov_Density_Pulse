#ifndef NORMA_H
#define NORMA_H

double l2_norma (double *u, double (*f) (double, double, double), int m, int row, double hx, double hy, double t, int offset = 0);
double c_norma (double *u, double (*f) (double, double, double), int m, int row, double t, double hx, double hy, int offset = 0);

void print_norms (int current_layer_t, double *density, double *pulse, int nodes_count, int half_nodes_count, int mx, double hx, double hy, double t);

#endif  //NORMA_H
