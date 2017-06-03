#include "stdio.h"
#include "scheme.h"
#include "solver.h"
#include "gas_params.h"
#include "time.h"

int main ()
{
    int it_max = 1;
    int it = 0;
    double time = 0.0;
    int mx = 10;
    int my = 10;
    int n = 10;
    gas_params params (1, 1, 1, mx, my, n);

    for (it = 0; it < it_max; it++)
    {
        time = clock ();

        scheme_Sokolov_Density_Pulse (&params);

        time = (clock() - time) / CLOCKS_PER_SEC;
        printf (">     Iter = %d     Time = %.4f\n", it, time);
        params.set_mutl_2 ();
    }
    return 0;
}
