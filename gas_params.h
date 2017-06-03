#ifndef GAS_PARAMS_H
#define GAS_PARAMS_H

struct gas_params
{
    double x;
    double y;
    double t;
    int mx;
    int my;
    int n;
    double h_x;
    double h_y;
    double tau;
    double mu;

    gas_params ();
    gas_params (double in_x, double in_y, double in_t, int in_mx, int in_my, int in_n);

    void set_mult_2 ();
    void update_steps ();
};

#endif // GAS_PARAMS_H
