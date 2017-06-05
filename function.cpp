#include "function.h"
#include "math.h"
#include "stdio.h"


double pulse_1 (double t, double x, double y)
{
    return sin (2 * M_PI * x) * sin (2 * M_PI * y) * exp (t);
}

double pulse_2 (double t, double x, double y)
{
    return sin (2 * M_PI * x) * sin (2 * M_PI * y) * exp (-t);
}

double density_ (double t, double x, double y)
{
    return (cos (2 * M_PI * x) + 3./2) * (sin (2 * M_PI * y) + 3./2) * exp (t);
}

double f1 (double t, double x, double y, double mu)
{
    double d = (3. + 2 * cos (2 * M_PI * x)) * (3. + 2 * sin (2 * M_PI * y)) * exp (2 * t);

    return (-52.6379 * mu * cos (2 * M_PI * x) * cos (2 * M_PI * y) + exp (2 * t) *
            (sin (2 * M_PI * x) * (-52.7788 - 35.1858 * sin (2 * M_PI * y)) *
             pow (exp (t) * (3./2 + cos (2 * M_PI * x)) * (3./2 + sin (2 * M_PI * y)), 0.4) -
             6.28319 * exp (2 * t) *
             (sin (2 * M_PI * x) - 6 * sin (4 * M_PI * x) - 3 * sin (6 * M_PI * x)) *
             sin (2 * M_PI * y) * sin (2 * M_PI * y) *
             (1.5 + sin (2 * M_PI * y)) +
             12 * exp (t) * sin (2 * M_PI * y) *
             (sin (4 * M_PI * x) * (0.5 + 1./3 * sin (2 * M_PI * y)) + sin (2 * M_PI * x) * (1.5 + sin (2 * M_PI * y))) +
             sin (2 * M_PI * x) * (368.465 * mu * sin (2 * M_PI * y) +
              (sin (4 * M_PI * x) * (18.8496 + 18.8496 * sin (2 * M_PI * y)) + sin (2 * M_PI * x) * (56.5487 + 56.5487 * sin (2 * M_PI * y))) * sin (4 * M_PI * y))
            )
           ) / d;
}

double f2 (double t, double x, double y, double mu)
{
    double d = (3. + 2 * cos (2 * M_PI * x)) * (3. + 2 * sin (2 * M_PI * y)) * exp (2 * t);

    return (exp (2 * t) * (52.7788 + 35.1858 * cos (2 * M_PI * x)) * cos (2 * M_PI * y) *
            pow (exp (t) * (3./2 + cos (2 * M_PI * x)) * (3./2 + sin (2 * M_PI * y)), 0.4) -
            52.6379 * exp (2 * t) * (mu * cos (2 * M_PI * x) * cos (2 * M_PI * y) +
                                     (sin (4 * M_PI * x) * (-1.0743 - 0.716197 * sin (2 * M_PI * y)) +
                                      sin (6 * M_PI * x) * (-0.537148 - 0.358099 * sin (2 * M_PI * y)) + sin (2 * M_PI * x) *
                                      (0.179049 + 0.119366 * sin (2 * M_PI * y))) * sin (2 * M_PI * y) * sin (2 * M_PI * y)) +
            sin (2 * M_PI * x) * (368.465 * mu * sin (2 * M_PI * y) + (sin (4 * M_PI * x) * (18.8496 + 18.8496 * sin (2 * M_PI * y)) + sin (2 * M_PI * x) *
                                                                       (56.5487 + 56.5487 * sin (2 * M_PI * y))) * sin (4 * M_PI * y))) / d;
}

void print_array (double *a, int m, int n)
{
    int i = 0, j = 0;
    printf ("<<<<<<<<<<<<\n");
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf (" %.4f ", a[i * n + j]);
        }
        printf ("\n");
    }
    printf ("<<<<<<<<<<<<\n");
}
