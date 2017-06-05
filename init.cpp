#include "init.h"
#include "function.h"

void init_density (double *d, gas_params *params, int half_nodes_count)
{
    int i = 0;
    int mx = params->mx;
    int my = params->my;
    double hx = params->h_x;
    double hy = params->h_y;
    double t = 0.0;

    for (i = 0; i < half_nodes_count; i++)
    {
        d[i] = density_ (t, i % (mx - 1) * hx + hx / 2, i / (mx - 1.) * hy + hy / 2);
    }
    print_array (d, my - 1, mx - 1);
}

void init_pulse (double *pulse, gas_params *params, int nodes_count)
{
    int i = 0;
    int mx = params->mx;
    int my = params->my;
    double hx = params->h_x;
    double hy = params->h_y;
    double *p1 = pulse;
    double *p2 = pulse + nodes_count;
    double t = 0.0;

    // y == 0
    for (i = 0; i < mx; i++)
    {
        p1[i] = 0.;
        p2[i] = 0.;
    }
    // 0 < y < maxY
    for (; i < nodes_count - mx; i++)
    {
        if (((i + 1) % mx == 0) || (i % mx == 0))
        {
            p1[i] = 0.;
            p2[i] = 0.;
            continue;
        }
        p1[i] = pulse_1(t, i % mx * hx, i / mx * hy);
        p2[i] = pulse_2(t, i % mx * hx, i / mx * hy);
    }
    // y == maxY
    for (; i < nodes_count; i++)
    {
        p1[i] = 0.;
        p2[i] = 0.;
    }
    print_array (p1, my, mx);
    print_array(p2, my, mx);
}
