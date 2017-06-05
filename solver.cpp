#include "solver.h"

#include "stdlib.h"
#include "memory.h"
#include "laspack/qmatrix.h"
#include "laspack/vector.h"
#include "laspack/rtc.h"
#include "laspack/itersolv.h"
#include "fill_in.h"

#define MAX_ITER 2000
#define EPS 1.e-8

static void get_result (double *result, Vector *x, int n)
{
    int i = 0;
    for (i = 0; i < n; i++)
        result[i] = V_GetCmp(x, i);
}

static void zero (double *pulse, int nodes_count, int mx)
{
    int i = 0;
    double *p1 = pulse;
    double *p2 = pulse + nodes_count;

    for (i = 0; i < mx; i++)
    {
        p1[i] = 0.0;
        p2[i] = 0.0;
    }
    for (; i < nodes_count - mx; i += mx)
    {
        p1[i] = 0.0;        p2[i] = 0.0;
        p1[i - 1] = 0.0;    p2[i - 1] = 0.0;
    }
    for (; i < nodes_count; i++)
    {
        p1[i] = 0.0;
        p2[i] = 0.0;
    }
}

void solve_density (double *density, double *old_density, int half_nodes_count, double *pulse, int nodes_count, gas_params *params, double t)
{
    QMatrix a;
    Vector x, b;

    Q_Constr (&a, "A", half_nodes_count, False, Rowws, Normal, True);
    V_Constr (&x, "x", half_nodes_count, Normal, True);
    V_Constr (&b, "b", half_nodes_count, Normal, True);

    SetRTCAccuracy(EPS);

    density_fill_matrix_and_rhs (&a, &b, &x, density, half_nodes_count, pulse, nodes_count, params, t);

    CGSIter(&a, &x, &b, MAX_ITER, SSORPrecond, 1);

    memcpy (old_density, density, sizeof (double) * half_nodes_count);
    get_result (density, &x, half_nodes_count);

    Q_Destr (&a);
    V_Destr (&x);
    V_Destr (&b);

    return;
}

void solve_pulse (double *density, double *old_density, double *pulse, double *old_pulse, int nodes_count, gas_params *params, double t)
{
    QMatrix a;
    Vector x, b;

    Q_Constr(&a, "A", 2 * nodes_count, False, Rowws, Normal, True);
    V_Constr(&x, "x", 2 * nodes_count, Normal, True);
    V_Constr(&b, "b", 2 * nodes_count, Normal, True);

    SetRTCAccuracy(EPS);

    pulse_fill_matrix_and_rhs (&a, &b, &x, density, old_density, pulse, nodes_count, params, t);

    CGSIter(&a, &x, &b, MAX_ITER, SSORPrecond, 1);

    memcpy (old_pulse, pulse, sizeof (double) * 2 * nodes_count);
    get_result (pulse, &x, nodes_count);
    zero (pulse, nodes_count, params->mx);

    Q_Destr (&a);
    V_Destr (&x);
    V_Destr (&b);

    return;
}
