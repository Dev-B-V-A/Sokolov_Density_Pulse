#include "scheme.h"
#include "solver.h"
#include "init.h"
#include "stdlib.h"
#include "report.h"
#include "stdio.h"
#include "norma.h"
#include "function.h"

void scheme_Sokolov_Density_Pulse (gas_params *params)
{
    const char *result_file = "report_scheme.tex";

    int nodes_count = params->mx * params->my;
    int half_nodes_count = (params->mx - 1) * (params->my - 1);
    int time_steps = params->n;
    int current_layer_t = 0;
    double *density;
    double *old_density;
    double *pulse;
    double *old_pulse;
    double mx = params->mx;
    double hx = params->h_x;
    double hy = params->h_y;
    double t = 0;
    double tau = params->tau;

    density = (double *)malloc (half_nodes_count * sizeof (double));
    old_density = (double *)malloc (half_nodes_count * sizeof (double));
    // pulse = [p1, p2]
    pulse = (double *)malloc (2 * nodes_count * sizeof (double));
    old_pulse = (double *)malloc (2 * nodes_count * sizeof (double));

    init_density (density, params, half_nodes_count);
    init_pulse (pulse, params, nodes_count);

    for (current_layer_t = 1; current_layer_t < time_steps; current_layer_t++)
    {
        t = current_layer_t * tau;
        printf (">      >       TimeStep = %d\n", current_layer_t);
        solve_density (density, old_density, half_nodes_count, pulse, nodes_count, params);

        solve_pulse (density, old_density, pulse, old_pulse, nodes_count, params, current_layer_t);

        printf ("Time_LAYER:    # %d\n", current_layer_t);
        printf ("Density          L2_norma = %.4f       C_norma = %.4f\n", l2_norma (density, density_, half_nodes_count, mx, hx, hy, t, 1),
                                                                           c_norma (density, density_, half_nodes_count, mx, t, hx, hy, 1));
        printf ("Pulse_1           L2_norma = %.4f       C_norma = %.4f\n", l2_norma (pulse, pulse_1, nodes_count, mx, hx, hy, t),
                                                                            c_norma (pulse, pulse_1, nodes_count, mx, t, hx, hy));
        printf ("Pulse_2           L2_norma = %.4f       C_norma = %.4f\n", l2_norma (pulse + nodes_count, pulse_2, nodes_count, mx, hx, hy, t),
                                                                            c_norma (pulse + nodes_count, pulse_2, nodes_count, mx, t, hx, hy));

    }

    record_report (result_file);

    free (density);
    free (pulse);
    free (old_density);
    free (old_pulse);
}
