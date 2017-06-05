#include "norma.h"
#include "math.h"
#include "stdio.h"
#include "function.h"

double c_norma (double *u, double (*f)(double, double, double), int m, int row, double t, double hx, double hy, int offset)
{
    double max = 0.0;
    double c = 0;
    double val = 0;
    int i = 0;
    for (i = 0; i < m; i++)
    {
        val = f (t, i % (row - offset) * hx + hx / 2 * offset, i / (row - offset) * hy + hy / 2);
        c = fabs (u[i] - val);
        if (c > max)
            max = c;
    }
    return max;
}

double l2_norma (double *u, double (*f)(double, double, double), int m, int row, double hx, double hy, double t, int offset)
{
    double sum = 0.0;
    double psum = 0.0;
    double val = 0;
    int i = 0;
    for (i = row; i < m - row; i++)
    {
        val = f (t, i % (row - offset) * hx + hx / 2 * offset, i / (row - offset) * hy + hy / 2);
        if (i % row == 0 || i % row == row - 1)
        {
            psum += u[i] * val;
            continue;
        }
        sum += u[i] * val;
    }
    return hx * hy * (sum + psum / 2);
}

void print_norms (int current_layer_t, double *density, double *pulse, int nodes_count, int half_nodes_count, int mx, double hx, double hy, double t)
{
    printf ("Time_LAYER:    # %d\n", current_layer_t);
    printf ("Density          L2_norma = %.4f       C_norma = %.4f\n", l2_norma (density, density_, half_nodes_count, mx, hx, hy, t, 1),
                                                                       c_norma (density, density_, half_nodes_count, mx, t, hx, hy, 1));
    printf ("Pulse_1           L2_norma = %.4f       C_norma = %.4f\n", l2_norma (pulse, pulse_1, nodes_count, mx, hx, hy, t),
                                                                        c_norma (pulse, pulse_1, nodes_count, mx, t, hx, hy));
    printf ("Pulse_2           L2_norma = %.4f       C_norma = %.4f\n", l2_norma (pulse + nodes_count, pulse_2, nodes_count, mx, hx, hy, t),
                                                                        c_norma (pulse + nodes_count, pulse_2, nodes_count, mx, t, hx, hy));
}
