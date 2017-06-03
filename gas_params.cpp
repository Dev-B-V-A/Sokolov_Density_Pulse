#include "gas_params.h"

#define MU 0

gas_params::gas_params()
{
    x = 0.0;
    y = 0.0;
    t = 0.0;
    mx = 0;
    my = 0;
    n = 0;
    h_x = 0.0;
    h_y = 0.0;
    tau = 0;
    mu = MU;
}

gas_params::gas_params (double in_x, double in_y, double in_t, int in_mx, int in_my, int in_n)
{
    x = in_x;
    y = in_y;
    t = in_t;
    mx = in_mx;
    my = in_my;
    n = in_n;
    mu = MU;
    update_steps();
}

void gas_params::set_mult_2 ()
{
    mx *= 2;
    my *= 2;
    n *= 2;
    update_steps();
}

void gas_params::update_steps ()
{
    h_x = x / (mx - 1);
    h_y = y / (my - 1);
    tau = t / (n - 1);
}
