#ifndef INIT_H
#define INIT_H

#include "gas_params.h"

void init_density (double *density, gas_params *params, int half_nodes_count);
void init_pulse (double *pulse, gas_params *params, int nodes_count);

#endif // INIT_H
