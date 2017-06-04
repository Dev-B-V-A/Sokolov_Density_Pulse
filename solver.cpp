#include "solver.h"

#include "stdlib.h"
#include "memory.h"
#include "qmatrix.h"
#include "vector.h"
#include "rtc.h"
#include "itersolv.h"
#include "fill_in.h"

#define MAX_ITER 2000
#define EPS 1.e-8

static void get_result (double *result, Vector *x, int n)
{
    int i = 0;
    for (i = 0; i < n; i++)
        result[i] = V_GetCmp(x, i);
}

void solve_density (double *density, double *old_density, int half_nodes_count, double *pulse, int nodes_count, gas_params *params)
{
    QMatrix a;
    Vector x, b;

    Q_Constr (&a, "A", half_nodes_count, False, Rowws, Normal, True);
    V_Constr (&x, "x", half_nodes_count, Normal, True);
    V_Constr (&b, "b", half_nodes_count, Normal, True);

    SetRTCAccuracy(EPS);

    density_fill_matrix_and_rhs (&a, &b, &x, density, half_nodes_count, pulse, nodes_count, params);

    CGSIter(&a, &x, &b, MAX_ITER, SSORPrecond, 1);

    get_result (density, &x, half_nodes_count);

    Q_Destr (&a);
    V_Destr (&x);
    V_Destr (&b);

    memcpy (old_density, density, sizeof (double) * half_nodes_count);
    return;
}

void solve_pulse (double *density, double *old_density, int half_nodes_count, double *pulse, double *old_pulse, int nodes_count, gas_params *params)
{
    QMatrix a;
    Vector x, b;

    Q_Constr(&a, "A", 2 * nodes_count, False, Rowws, Normal, True);
    V_Constr(&x, "x", 2 * nodes_count, Normal, True);
    V_Constr(&b, "b", 2 * nodes_count, Normal, True);

    SetRTCAccuracy(EPS);

    pulse_fill_matrix_and_rhs (&a, &b, &x, density, old_density, half_nodes_count, pulse, nodes_count, params);

    CGSIter(&a, &x, &b, MAX_ITER, SSORPrecond, 1);

    get_result (pulse, &x, nodes_count);

    Q_Destr (&a);
    V_Destr (&x);
    V_Destr (&b);

    memcpy (old_pulse, pulse, sizeof (double) * 2 * nodes_count);
    return;
}
