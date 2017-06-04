#include "solver.h"
#include "stdlib.h"
#include "memory.h"
#include "qmatrix.h"
#include "rtc.h"
#include "itersolv.h"

#define MAX_ITER 2000

void solve_density (double *density, double *old_density, int half_nodes_count, double *pulse, double *old_pulse, int nodes_count, gas_params *params)
{
    QMatrix a;
    Vector x, b;

    Q_Constr(&a, "A", half_nodes_count, false, Rowws, Normal, true);
    V_Constr(&x, "x", half_nodes_count, Normal, true);
    V_Constr(&b, "b", half_nodes_count, Normal, true);

    fill_matrix (&a, density, half_nodes_count, pulse, nodes_count, params);

    fill_init_x (&x, density, half_nodes_count, pulse, nodes_count, params);

    fill_rhs (&b, density, half_nodes_count, pulse, nodes_count, params);

    CGSIter(&a, &x, &b, MAX_ITER, SSORPrecond, 1);

    memcpy (old_density, density, sizeof (double) * half_nodes_count);
    return;
}

void solve_pulse (double *density, double *old_density, int half_nodes_count, double *pulse, double *old_pulse, int nodes_count, gas_params *params)
{

    memcpy (old_pulse, pulse, sizeof (double) * 2 * nodes_count);
    return;
}
