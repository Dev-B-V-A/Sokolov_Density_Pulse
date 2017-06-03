#include "init.h"

void init_density (double *density, gas_params *params, int half_nodes_count)
{
    int i = 0;
    int j = 0;
    int mx = params->mx;
    int my = params->my;
    double hx = params->h_x;
    double hy = params->h_y;
    for (i = 0; i < my; i++)
    {
        if (i == 0)
        {
            for (j = 0; j < mx; j++)
            {

            }
        }
        if (i == my - 1)
        {

        }
    }
}

void init_pulse (double *pulse, gas_params *params, int nodes_count)
{

}
