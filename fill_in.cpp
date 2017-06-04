#include "fill_in.h"

#include "qmatrix.h"
#include "itersolv.h"

static double p1_tilde (double *p1, int current_node, int mx)
{
    return (p1[current_node] + p1[current_node + mx]) / 2.;
}

static double p2_tilde (double *p2, int current_node)
{
    return (p2[current_node] + p2[current_node + 1]) / 2.;
}

//static double h_tilde (double *h, int current_node, int mx)
//{
//    return (h[current_node] + h[current_node - mx] + h[current_node - 1] + h[current_node - mx - 1]) / 4.;
//}

//static double h1_tilde (double *h, int current_node, int mx)
//{
//    return (h[current_node] + h[current_node - mx]) / 2.;
//}

//static double h2_tilde (double *h, int current_node)
//{
//    return (h[current_node] + h[current_node - 1]) / 2.;
//}

static double ah0 (int current_node, double *pulse, double hx, double hy, double tau, int nodes_count, int mx)
{
    double *p1 = pulse;
    double *p2 = pulse + nodes_count;
    return 1./tau +
           (  1./2/hx * (p1_tilde(p1, current_node + 1, mx) +
                          fabs (p1_tilde (p1, current_node + 1, mx)) -
                          p1_tilde (p1, current_node, mx) +
                          fabs (p1_tilde (p1, current_node, mx))) +
              1./2/hy * (p2_tilde (p2, current_node + mx) +
                         fabs (p2_tilde (p2, current_node + mx)) -
                         p2_tilde (p2, current_node) +
                         fabs (p2_tilde (p2, current_node)))
           );
}

static double ahr (int current_node, double *pulse, double hx, int mx)
{
    double *p1 = pulse;
    return 1./2/hx *
            (p1_tilde (p1, current_node + 1, mx) -
             fabs (p1_tilde (p1, current_node, mx))
            );
}

static double ahup (int current_node, double *pulse, int nodes_count, double hy, int mx)
{
    double *p2 = pulse + nodes_count;
    return 1./2/hy *
            (p2_tilde (p2, current_node + mx) -
             fabs (p2_tilde (p2, current_node + mx))
            );

}

static double ahl (int current_node, double *pulse, double hx, int mx)
{
    double *p1 = pulse;
    return -1./2/hx *
            (p1_tilde (p1, current_node, mx) +
             fabs (p1_tilde (p1, current_node, mx))
             );
}

static double ahd (int current_node, double *pulse, double hy)
{
    double *p2 = pulse;
    return -1./2/hy *
            (p2_tilde (p2, current_node) +
             fabs (p2_tilde (p2, current_node))
            );
}

static int get_up (int current_node, int mx)
{
    return current_node + mx;
}

static int get_right (int current_node)
{
    return current_node + 1;
}

static int get_left (int current_node)
{
    return current_node - 1;
}

static int get_down (int current_node, int mx)
{
    return current_node - mx;
}

void density_fill_matrix_and_rhs (QMatrix *a, Vector *b, Vector *x, double *density, int half_nodes_count, double *pulse, int nodes_count, gas_params *params)
{
    int current_node = 0;
    double tmp = 0;
    int nz_angle = 3;
    int nz_side = 4;
    int nz = 5;
    int curr_nz = 0;
    double tau = params->tau;
    int mx = params->mx;
    double hx = params->h_x;
    double hy = params->h_y;

    // Left Down Angle
    {
        Q_SetLen (a, current_node, nz_angle);
        // 00
        tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
        curr_nz = 0;
        Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
        // 0R
        tmp = ahr (current_node, pulse, hx, mx);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_right(current_node), tmp);
        // 0U
        tmp = ahup (current_node, pulse, nodes_count, hy, mx);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_up(current_node, mx), tmp);

        V_SetCmp (b, current_node, density[current_node] / tau);
        V_SetCmp (x, current_node, current_node % 2);

        current_node++;
    }

    // Down side
    {
        for (; current_node < mx - 1; current_node++)
        {
            Q_SetLen (a, current_node, nz_side);
            // 00
            tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
            curr_nz = 0;
            Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
            // 0R
            tmp = ahr (current_node, pulse, hx, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_right(current_node), tmp);
            // 0U
            tmp = ahup (current_node, pulse, nodes_count, hy, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_up(current_node, mx), tmp);
            // 0L
            tmp = ahl (current_node, pulse, hx, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

            V_SetCmp (b, current_node, density[current_node] / tau);
            V_SetCmp (x, current_node, current_node % 2);
        }

    }

    // Right Down Angle
    {
        Q_SetLen (a, current_node, nz_angle);
        // 00
        tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
        curr_nz = 0;
        Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
        // 0L
        tmp = ahl (current_node, pulse, hx, mx);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);
        // 0U
        tmp = ahup (current_node, pulse, nodes_count, hy, mx);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_up(current_node, mx), tmp);

        V_SetCmp (b, current_node, density[current_node] / tau);
        V_SetCmp (x, current_node, current_node % 2);

        current_node++;
    }

    // Input points
    {
        for (; current_node < half_nodes_count - mx + 1; current_node++)
        {
            // Left side
            if (current_node % (mx - 1) == 0)
            {
                Q_SetLen (a, current_node, nz_side);
                // 00
                tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
                curr_nz = 0;
                Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
                // 0R
                tmp = ahr (current_node, pulse, hx, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_right(current_node), tmp);
                // 0U
                tmp = ahup (current_node, pulse, nodes_count, hy, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_up(current_node, mx), tmp);
                // 0D
                tmp = ahd (current_node, pulse, hy);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

                V_SetCmp (b, current_node, density[current_node] / tau);
                V_SetCmp (x, current_node, current_node % 2);
                continue;
            }

            // Right side
            if (current_node % (mx - 1) == mx - 2)
            {
                Q_SetLen (a, current_node, nz_side);
                // 00
                tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
                curr_nz = 0;
                Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
                // 0D
                tmp = ahd (current_node, pulse, hy);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);
                // 0U
                tmp = ahup (current_node, pulse, nodes_count, hy, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_up(current_node, mx), tmp);
                // 0L
                tmp = ahl (current_node, pulse, hx, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

                V_SetCmp (b, current_node, density[current_node] / tau);
                V_SetCmp (x, current_node, current_node % 2);

                continue;
            }

            Q_SetLen (a, current_node, nz);
            // 00
            tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
            curr_nz = 0;
            Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
            // 0R
            tmp = ahr (current_node, pulse, hx, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_right(current_node), tmp);
            // 0U
            tmp = ahup (current_node, pulse, nodes_count, hy, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_up(current_node, mx), tmp);
            // 0L
            tmp = ahl (current_node, pulse, hx, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);
            // 0D
            tmp = ahd (current_node, pulse, hy);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

            V_SetCmp (b, current_node, density[current_node] / tau);
            V_SetCmp (x, current_node, current_node % 2);
        }

    }

    // Left Up Angle
    {
        Q_SetLen (a, current_node, nz_angle);
        // 00
        tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
        curr_nz = 0;
        Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
        // 0R
        tmp = ahr (current_node, pulse, hx, mx);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_right(current_node), tmp);
        // 0D
        tmp = ahd (current_node, pulse, hy);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

        V_SetCmp (b, current_node, density[current_node] / tau);
        V_SetCmp (x, current_node, current_node % 2);

        current_node++;
    }

    // Up side
    {
        for (; current_node < mx - 1; current_node++)
        {
            Q_SetLen (a, current_node, nz_side);
            // 00
            tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
            curr_nz = 0;
            Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
            // 0R
            tmp = ahr (current_node, pulse, hx, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_right(current_node), tmp);
            // 0D
            tmp = ahd (current_node, pulse, hy);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);
            // 0L
            tmp = ahl (current_node, pulse, hx, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

            V_SetCmp (b, current_node, density[current_node] / tau);
            V_SetCmp (x, current_node, current_node % 2);
        }

    }

    // Up Right Angle
    {
        Q_SetLen (a, current_node, nz_angle);
        // 00
        tmp = ah0 (current_node, pulse, hx, hy, tau, nodes_count, mx);
        curr_nz = 0;
        Q_SetEntry (a, current_node, curr_nz, current_node, tmp);
        // 0L
        tmp = ahl (current_node, pulse, hx, mx);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);
        // 0D
        tmp = ahd (current_node, pulse, hy);
        curr_nz++;
        Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

        V_SetCmp (b, current_node, density[current_node] / tau);
        V_SetCmp (x, current_node, current_node % 2);

        current_node++;
    }
}

void pulse_fill_matrix_and_rhs (QMatrix */*a*/, Vector */*b*/, Vector */*x*/, double */*density*/, double */*old_density*/,
                                int /*half_nodes_count*/, double */*pulse*/, int /*nodes_count*/, gas_params */*params*/)
{
    // Left Down Angle

    // Down side

    // Right Down Angle

    // Input points

    // Left side

    // Right side

    // Left Up Angle

    // Up side

    // Up Right Angle
}
