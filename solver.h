#ifndef SOLVER_H
#define SOLVER_H

#include "gas_params.h"

void solve_density (double *density, double *old_density, int half_nodes_count,
                    double *pulse, int nodes_count, gas_params *params, double t);

void solve_pulse (double *density, double *old_density,
                  double *pulse, double *old_pulse, int nodes_count, gas_params *params, double t);

#endif // SOLVER_H
