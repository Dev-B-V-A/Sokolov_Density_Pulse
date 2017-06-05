#include "norma.h"
#include "math.h"

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
