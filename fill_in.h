#ifndef FILL_IN_H
#define FILL_IN_H

#include "laspack/qmatrix.h"
#include "gas_params.h"

void density_fill_matrix_and_rhs (QMatrix *a, Vector *b, Vector *x, double *density, int half_nodes_count,
                                  double *pulse, int nodes_count, gas_params *params);

void pulse_fill_matrix_and_rhs (QMatrix *a, Vector *b, Vector *x, double *density, double *old_density, double *pulse, int nodes_count, gas_params *params, double t);

#endif // FILL_IN_H
