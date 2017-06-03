#include "scheme.h"
#include "solver.h"

void scheme_Sokolov_Density_Pulse (gas_params *params)
{
    const char *result_file = "report_scheme.tex";

    int nodes_count = params->mx * params->my;
    int half_nodes_count = (params->mx + 1) * (params->my + 1);
    int time_steps = params->n;
    int current_layer_t = 0;
    double *density;
    double *old_density;
    double *pulse;
    double *old_pulse;

    density = (double *)malloc (half_nodes_count * sizeof (double));
    old_density = (double *)malloc (half_nodes_count * sizeof (double));
    pulse = (double *)malloc (2 * nodes_count * sizeof (double));
    old_pulse = (double *)malloc (2 * nodes_count * sizeof (double));

    init_H (density, params, half_nodes_count);
    init_V (pulse, params, nodes_count * 2);

    for (current_layer_t = 0; current_layer_t < time_steps; current_layer_t++)
    {
        solve_density ();

        solve_pulse();
    }

    record_report (result_file);

    free (density);
    free (pulse);
    free (old_density);
    free (old_pulse);
}
