#include "fill_in.h"

#include "laspack/qmatrix.h"
#include "laspack/itersolv.h"

#include "function.h"

#define GAMMA 1.4

static double p1_tilde (double *p1, int current_node, int mx)
{
    return (p1[current_node] + p1[current_node + mx]) / 2.;
}

static double p2_tilde (double *p2, int current_node)
{
    return (p2[current_node] + p2[current_node + 1]) / 2.;
}

//static double get_h_from_p (double *h, int current_node, int mx, int my)
//{
//    int i = current_node / mx;
//    int j = current_node % mx;

//    if (j == (mx - 1) || i == my - 1)
//        return 0;
//    return h[i * (mx - 1) + j];
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

void density_fill_matrix_and_rhs (QMatrix *a, Vector *b, Vector *x, double *density, int half_nodes_count,
                                  double *pulse, int nodes_count, gas_params *params)
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

static double h_tilde (double *h, int current_node, int mx)
{
    if (current_node / mx == 0 && current_node % mx == 0)
        return 1./4 * h[current_node];
    if (current_node / mx == 0)
        return 1./4 * (h[current_node] + h[current_node - 1]);
    if (current_node % mx == 0)
        return 1./4 * (h[current_node] + h[current_node - mx]);

    return (h[current_node] + h[current_node - mx] + h[current_node - 1] + h[current_node - mx - 1]) / 4.;
}

static double h1_tilde (double *h, int current_node, int mx)
{
    if (current_node / mx == 0)
        return h[current_node] / 2;
    return (h[current_node] + h[current_node - mx]) / 2.;
}

static double h2_tilde (double *h, int current_node, int mx)
{
    if (current_node % mx == 0)
        return h[current_node] / 2;
    return (h[current_node] + h[current_node - 1]) / 2.;
}

static double get_p_left (double *p1, int current_node, int mx)
{
    if (current_node % mx == 0)
        return 0;

    return p1[current_node - 1];
}

static double get_h1_left_2 (double *density, int current_node, int mx)
{
    if (current_node % (mx - 1) == 0 || current_node % (mx - 1) == 1)
        return 0;
    return h1_tilde (density, current_node - 2, mx);
}

static double get_h1_left (double *density, int current_node, int mx)
{
    if (current_node % (mx - 1) == 0)
        return 0;
    return h1_tilde (density, current_node - 1, mx);
}

static double get_h1_right (double *density, int current_node, int mx)
{
    if (current_node % (mx - 1) == (mx - 2))
        return 0;
    return h1_tilde (density, current_node + 1, mx);
}

static double get_p_right (double *p1, int current_node, int mx)
{
    if (current_node % mx == (mx - 1))
        return 0;
    return p1[current_node + 1];
}

static double get_p_down (double *p1, int current_node, int mx)
{
    if (current_node / mx == 0)
        return 0;
    return p1[current_node - mx];
}

static double get_h2_down_2 (double *density, int current_node, int mx)
{
    if (current_node / mx == 0 || current_node / mx == 1)
        return 0;
    return h2_tilde (density, current_node - 2 * mx, mx);
}

static double get_h2_down (double *density, int current_node, int mx)
{
    if (current_node / mx == 0)
        return 0;
    return h2_tilde (density, current_node - mx, mx);
}

static double get_p_up (double *pulse, int current_node, int mx, int my)
{
    if (current_node / mx == my - 1)
        return 0;
    return pulse[current_node + mx];
}

static double get_h2_up (double *density, int current_node, int mx, int my)
{
    if (current_node / mx == my - 1)
        return 0;
    return h2_tilde (density, current_node + mx, mx);
}

static double ap10 (int current_node,  double *density, double *pulse,
                    double hx, double hy, double tau, double mu, int mx)
{
    double *p1 = pulse;
    return h_tilde (density, current_node, mx) / tau +
            2 * mu * (4./3 / hx / hx + 1./ hy / hy) +
            1./4/hx * (
                (fabs (get_p_left (p1, current_node, mx) - get_p_left(p1, current_node, mx) +
                 fabs (p1[current_node]) + p1[current_node]) *
                 get_h1_left(density, current_node, mx) +
                 (fabs (get_p_right (p1, current_node, mx)) + get_p_right (p1, current_node, mx) +
                  fabs (p1[current_node]) - p1[current_node]) *
                 h1_tilde (density, current_node, mx)));
}

static double ap1l (int current_node, double *density, double *pulse, double hx,  double mu, int mx)
{
    double *p1 = pulse;
    return -(1./4/hx *
              ((fabs (get_p_left(p1, current_node, mx)) + get_p_left(p1, current_node, mx)) * get_h1_left_2 (density, current_node, mx) +
               (fabs (p1[current_node]) + p1[current_node]) * get_h1_left(density, current_node, mx)
             ) +
             mu * 4. / 3 / hx / hx);
}

static double ap1r (int current_node, double *density, double *pulse, double hx, double mu, int mx)
{
    double *p1 = pulse;
    return -(1./4/hx *
             ((fabs (p1[current_node]) - p1[current_node]) * h1_tilde (density, current_node, mx) +
              (fabs (get_p_right (p1, current_node, mx)) - get_p_right(p1, current_node, mx)) * get_h1_right (density, current_node, mx) +
              mu * 4./3/hx/hx
                 ));
}

static double ap2d (int current_node, double *density, double *pulse, double hy, int mx)
{
    double *p1 = pulse;
    return -(1./4/hy *
             ((fabs (get_p_down (p1, current_node, mx)) + get_p_down (p1, current_node, mx)) *
              get_h2_down_2 (density, current_node, mx) +
              (fabs (p1[current_node]) + p1[current_node]) *
              get_h2_down (density, current_node, mx)
             )
            );
}

static double ap20 (int current_node, double *density, double *pulse, double hy, int mx, int my)
{
    double *p1 = pulse;
    return 1./4/hy *
            ((fabs (get_p_down (p1, current_node, mx)) - get_p_down(p1, current_node, mx) + fabs (p1[current_node]) + p1[current_node]) *
             get_h2_down (density, current_node, mx) +
             (fabs (get_p_up (p1, current_node, mx, my)) + get_p_up (p1, current_node, mx, my) + fabs (p1[current_node]) - p1[current_node]) *
             h2_tilde (density, current_node, mx)
                );
}

static double ap2u (int current_node, double *density, double *pulse, double hy, int mx, int my)
{
    double *p1 = pulse;
    return -(1./4/hy *
             ((fabs (p1[current_node]) - p1[current_node]) *
              h2_tilde (density, current_node, mx) +
              (fabs (get_p_up (p1, current_node, mx, my)) - get_p_up (p1, current_node, mx ,my)) *
              get_h2_up (density, current_node, mx, my))
            );
}

static double ap1d (double mu, double hy)
{
    return -mu /hy/hy;
}

static double ap1u (double mu, double hy)
{
    return -mu/hy/hy;
}

static double get_p_ld (double *p, int current_node, int mx)
{
    if (current_node / mx == 0 || current_node % mx == 0)
        return 0;
    return p[current_node - mx - 1];
}

static double get_p_lu (double *p, int current_node, int mx, int my)
{
    if (current_node / mx == my - 1 || current_node % mx == 0)
        return 0;
    return p[current_node + mx - 1];
}

static double get_p_rd (double *p, int current_node, int mx)
{
    if (current_node / mx == 0 || current_node % mx == mx - 1)
        return 0;
    return p[current_node - mx + 1];
}

static double get_p_ru (double *p, int current_node, int mx, int my)
{
    if (current_node / mx == my - 1 || current_node % mx == mx - 1)
        return 0;
    return p[current_node + mx + 1];
}

static double a_rhs (int current_node, double *density, double *old_density,
                     double *pulse, double tau, double gamma, double mu, double hx, double hy,
                     int mx, int my, int nodes_count, double t)
{
    double *p1 = pulse;
    double *p2 = pulse + nodes_count;
    double x = current_node % mx * hx;
    double y = current_node / mx * hy;
    return h_tilde (old_density, current_node, mx) * p1[current_node] / tau +
            mu / 12 /hx/hy *
            (get_p_ld (p2, current_node, mx) - get_p_lu (p2, current_node, mx,my) -
             get_p_rd (p2, current_node, mx) + get_p_ru (p2, current_node, mx, my)) +
            f1 (t, x, y, mu) * h_tilde (density, current_node, mx) -
            gamma / (gamma - 1) * h_tilde (density, current_node, mx) / hx *
            (pow (h1_tilde (density, current_node, mx), gamma - 1) - pow (get_h1_left (density, current_node, mx), gamma - 1));
}

static double bp20 (int current_node, double *p2, double *density, int mx, int my, double hx, double hy, double tau, double mu)
{
    return h_tilde (density, current_node, mx) / tau +
            1./4/hy *
            ((fabs (get_p_down (p2, current_node, mx)) - get_p_down (p2, current_node, mx) + fabs (p2[current_node]) + p2[current_node]) *
             get_h2_down(density, current_node, mx) +
             (fabs (get_p_up (p2, current_node, mx, my)) + get_p_up (p2, current_node, mx, my) + fabs (p2[current_node]) - p2[current_node]) *
             h2_tilde (density, current_node, mx)
            ) +
            2 * mu * (1/hx/hx + 4./3/hy/hy);

}

static double bp1l (int current_node, double *p2, double *density, int mx, int my, double hx, double mu)
{
    return -1./4/hx *
             ((fabs (get_p_left (p2, current_node, mx)) + get_p_left (p2, current_node, mx)) *
              get_h1_left_2(density, current_node, mx) +
              (fabs (p2[current_node]) + p2[current_node]) *
              get_h1_left (density, current_node, mx));
}

static double bp10 (int current_node, double *p2, double *density, int mx, int my, double hx)
{
    return 1./4/hx *
            ((fabs (get_p_left (p2, current_node, mx)) - get_p_left (p2, current_node, mx) + fabs (p2[current_node]) + p2[current_node]) *
             get_h1_left(density, current_node, mx) +
             (fabs (get_p_right (p2, current_node, mx)) + get_p_right(p2, current_node, mx) + fabs (p2[current_node]) - p2[current_node]) *
             h1_tilde (density, current_node, mx));
}

static double bp1r (int current_node, double *p2, double *density, int mx, int my, double hx, double mu)
{
    return -1./4/hx *
            ((fabs (p2[current_node]) - p2[current_node]) *
             h1_tilde (density, current_node, mx) +
             (fabs (get_p_right (p2, current_node, mx)) - get_p_right (p2, current_node, mx)) *
             get_h1_right(density, current_node, mx));
}

static double bp2d (int current_node, double *p2, double *density, int mx, int my, double hy, double mu)
{
    return -1./4/hy *
            ((fabs (get_p_down (p2, current_node, mx)) + get_p_down(p2, current_node, mx)) *
             get_h2_down_2 (density, current_node, mx) +
             (fabs (p2[current_node]) + p2[current_node]) *
             get_h2_down(density, current_node, mx)) -
            mu* 4./3 /hy/hy;
}
static double bp2u (int current_node, double *p2, double *density, int mx, int my, double hy, double mu)
{
    return -1./4/hy *
            ((fabs (p2[current_node]) - p2[current_node]) *
             h2_tilde (density, current_node, mx) +
             (fabs (get_p_up (p2, current_node, mx, my)) - get_p_up (p2, current_node, mx, my)) *
             get_h2_up(density, current_node, mx, my)) -
            mu * 4./3 /hy/hy;
}

static double bp2l (double mu, double hx)
{
    return mu /hx/hx;
}

static double bp2r (double mu, double hx)
{
    return mu/hx/hx;
}

static double b_rhs (int current_node, double *pulse, double *old_density, double *density,
                     double mu, double hx, double hy, int mx, int my, int nodes_count, double tau, double gamma, double t)
{
    double *p1 = pulse;
    double *p2 = pulse + nodes_count;
    double x = current_node % mx * hx;
    double y = current_node / mx * hy;
    return h_tilde (old_density, current_node, mx) * p2[current_node] / tau -
            gamma / (gamma - 1) / hy * h_tilde(density, current_node, mx) *
            (pow (h2_tilde (density, current_node, mx), gamma - 1) - pow (get_h2_down (density, current_node, mx), gamma - 1)) +
             mu /12/hx/hy * (get_p_ld(p1, current_node, mx) - get_p_lu (p1, current_node, mx, my) -
                              get_p_rd (p1, current_node, mx) + get_p_ru (p1, current_node, mx, my)) +
             f2 (t, x, y, mu) * h_tilde(density, current_node, mx);
}

void pulse_fill_matrix_and_rhs (QMatrix *a, Vector *b, Vector *x, double *density, double *old_density,
                                double *pulse, int nodes_count, gas_params *params, double t)
{
    int current_node = 0;
    int nz = 8;
    int nz_angle = 5;
    int nz_side_vert1, nz_side_hor2 = 6;
    int nz_side_hor1, nz_side_vert2 = 7;
    int curr_nz= 0;
    double tmp = 0;
    double mu = params->mu;
    double tau = params->tau;
    double *p2 = pulse + nodes_count;

    int mx = params->mx;
    int my = params->my;
    double hx = params->h_x;
    double hy = params->h_y;
    // Left Down Angle
    {
        Q_SetLen (a, current_node, nz_angle);
        {
            tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
            curr_nz = 0;
            Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

            tmp = ap1r(current_node, density, pulse, hx, mu, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_right (current_node), tmp);

            tmp = ap1u (mu, hy);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_up (current_node, mx), tmp);

            tmp = ap20 (current_node, density, pulse, hy, mx ,my);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

            tmp = ap2u(current_node, density, pulse, hy, mx, my);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, nodes_count + get_up (current_node, mx), tmp);

            V_SetCmp(x, current_node, current_node % 2);
            tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
            V_SetCmp(b, current_node, tmp);
        }

        Q_SetLen(a, current_node + nodes_count, nz_angle);
        {
            tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
            curr_nz = 0;
            Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

            tmp = bp2u(current_node, p2, density, mx, my, hy, mu);
            curr_nz++;
            Q_SetEntry(a, current_node + nodes_count, curr_nz, nodes_count + get_up (current_node, mx), tmp);

            tmp = bp2r(mu, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_right (current_node), tmp);

            tmp = bp10(current_node, p2, density, mx, my, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

            tmp = bp1r(current_node, p2, density, mx, my, hx, mu);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, get_right (current_node), tmp);

            V_SetCmp(x, current_node + nodes_count, current_node % 2);
            tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
            V_SetCmp(b, nodes_count + current_node, tmp);
        }
        current_node++;
    }

    // Down side
    {
        for (; current_node < mx; current_node++)
        {
            Q_SetLen(a, current_node, nz_side_vert1);
            {
                tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
                curr_nz = 0;
                Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

                tmp = ap1r(current_node, density, pulse, hx, mu, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_right (current_node), tmp);

                tmp = ap1u (mu, hy);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_up (current_node, mx), tmp);

                tmp = ap1l(current_node, density, pulse, hx, mu, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

                tmp = ap20 (current_node, density, pulse, hy, mx ,my);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

                tmp = ap2u(current_node, density, pulse, hy, mx, my);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                V_SetCmp(x, current_node, current_node % 2);
                tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
                V_SetCmp(b, current_node, tmp);

            }
            Q_SetLen(a, nodes_count + current_node, nz_side_vert2);
            {
                tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
                curr_nz = 0;
                Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

                tmp = bp2u(current_node, p2, density, mx, my, hy, mu);
                curr_nz++;
                Q_SetEntry(a, current_node + nodes_count, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                tmp = bp2r(mu, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_right (current_node), tmp);

                tmp = bp2l (mu, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_left (current_node), tmp);

                tmp = bp10(current_node, p2, density, mx, my, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

                tmp = bp1r(current_node, p2, density, mx, my, hx, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, get_right (current_node), tmp);

                tmp = bp1l(current_node, p2, density, mx, my, hx, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, get_left (current_node), tmp);

                V_SetCmp(x, current_node + nodes_count, current_node % 2);
                tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
                V_SetCmp(b, nodes_count + current_node, tmp);
            }
        }
    }

    // Right Down Angle
    {
        Q_SetLen (a, current_node, nz_angle);
        {
            tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
            curr_nz = 0;
            Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

            tmp = ap1l(current_node, density, pulse, hx, mu, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

            tmp = ap1u (mu, hy);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_up (current_node, mx), tmp);

            tmp = ap20 (current_node, density, pulse, hy, mx ,my);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

            tmp = ap2u(current_node, density, pulse, hy, mx, my);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, nodes_count + get_up (current_node, mx), tmp);

            V_SetCmp(x, current_node, current_node % 2);
            tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
            V_SetCmp(b, current_node, tmp);
        }

        Q_SetLen(a, current_node + nodes_count, nz_angle);
        {
            tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
            curr_nz = 0;
            Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

            tmp = bp2u(current_node, p2, density, mx, my, hy, mu);
            curr_nz++;
            Q_SetEntry(a, current_node + nodes_count, curr_nz, nodes_count + get_up (current_node, mx), tmp);

            tmp = bp2l (mu, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_left (current_node), tmp);

            tmp = bp10(current_node, p2, density, mx, my, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

            tmp = bp1l(current_node, p2, density, mx, my, hx, mu);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, get_left (current_node), tmp);

            V_SetCmp(x, current_node + nodes_count, current_node % 2);
            tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
            V_SetCmp(b, nodes_count + current_node, tmp);
        }
        current_node++;
    }

    // Input points
    {
        for (; current_node < nodes_count - mx; current_node++)
        {
            // Left side
            if (current_node % mx == 0)
            {
                Q_SetLen(a, current_node, nz_side_hor1);
                {
                    tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
                    curr_nz = 0;
                    Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

                    tmp = ap1r(current_node, density, pulse, hx, mu, mx);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, get_right (current_node), tmp);

                    tmp = ap1u (mu, hy);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, get_up (current_node, mx), tmp);

                    tmp = ap1d (mu, hy);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

                    tmp = ap20 (current_node, density, pulse, hy, mx ,my);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

                    tmp = ap2u(current_node, density, pulse, hy, mx, my);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                    tmp = ap2d(current_node, density, pulse, hy, mx);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                    V_SetCmp(x, current_node, current_node % 2);
                    tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
                    V_SetCmp(b, current_node, tmp);

                }
                Q_SetLen(a, nodes_count + current_node, nz_side_hor2);
                {
                    tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
                    curr_nz = 0;
                    Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

                    tmp = bp2u(current_node, p2, density, mx, my, hy, mu);
                    curr_nz++;
                    Q_SetEntry(a, current_node + nodes_count, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                    tmp = bp2r(mu, hx);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_right (current_node), tmp);

                    tmp = bp2d (current_node, p2, density, mx, my, hy, mu);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                    tmp = bp10(current_node, p2, density, mx, my, hx);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

                    tmp = bp1r(current_node, p2, density, mx, my, hx, mu);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, get_right (current_node), tmp);

                    V_SetCmp(x, current_node + nodes_count, current_node % 2);
                    tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
                    V_SetCmp(b, nodes_count + current_node, tmp);
                }
                continue;
            }

            // Right side
            if (current_node % mx == mx - 1)
            {
                Q_SetLen(a, current_node, nz_side_hor1);
                {
                    tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
                    curr_nz = 0;
                    Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

                    tmp = ap1l(current_node, density, pulse, hx, mu, mx);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

                    tmp = ap1u (mu, hy);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, get_up (current_node, mx), tmp);

                    tmp = ap1d (mu, hy);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

                    tmp = ap20 (current_node, density, pulse, hy, mx ,my);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

                    tmp = ap2u(current_node, density, pulse, hy, mx, my);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                    tmp = ap2d(current_node, density, pulse, hy, mx);
                    curr_nz++;
                    Q_SetEntry(a, current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                    V_SetCmp(x, current_node, current_node % 2);
                    tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
                    V_SetCmp(b, current_node, tmp);

                }
                Q_SetLen(a, nodes_count + current_node, nz_side_hor2);
                {
                    tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
                    curr_nz = 0;
                    Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

                    tmp = bp2u(current_node, p2, density, mx, my, hy, mu);
                    curr_nz++;
                    Q_SetEntry(a, current_node + nodes_count, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                    tmp = bp2l(mu, hx);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_right (current_node), tmp);

                    tmp = bp2d (current_node, p2, density, mx, my, hy, mu);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                    tmp = bp10(current_node, p2, density, mx, my, hx);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

                    tmp = bp1l(current_node, p2, density, mx, my, hx, mu);
                    curr_nz++;
                    Q_SetEntry(a, nodes_count + current_node, curr_nz, get_right (current_node), tmp);

                    V_SetCmp(x, current_node + nodes_count, current_node % 2);
                    tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
                    V_SetCmp(b, nodes_count + current_node, tmp);
                }
                continue;
            }
            Q_SetLen(a, current_node, nz);
            {
                tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
                curr_nz = 0;
                Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

                tmp = ap1r(current_node, density, pulse, hx, mu, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_right (current_node), tmp);

                tmp = ap1u (mu, hy);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_up (current_node, mx), tmp);

                tmp = ap1l(current_node, density, pulse, hx, mu, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

                tmp = ap1d (mu, hy);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

                tmp = ap20 (current_node, density, pulse, hy, mx ,my);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

                tmp = ap2u(current_node, density, pulse, hy, mx, my);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                tmp = ap2d(current_node, density, pulse, hy, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                V_SetCmp(x, current_node, current_node % 2);
                tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
                V_SetCmp(b, current_node, tmp);

            }
            Q_SetLen(a, nodes_count + current_node, nz);
            {
                tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
                curr_nz = 0;
                Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

                tmp = bp2u(current_node, p2, density, mx, my, hy, mu);
                curr_nz++;
                Q_SetEntry(a, current_node + nodes_count, curr_nz, nodes_count + get_up (current_node, mx), tmp);

                tmp = bp2r(mu, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_right (current_node), tmp);

                tmp = bp2l (mu, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_left (current_node), tmp);

                tmp = bp2d (current_node, p2, density, mx, my, hy, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                tmp = bp10(current_node, p2, density, mx, my, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

                tmp = bp1r(current_node, p2, density, mx, my, hx, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, get_right (current_node), tmp);

                tmp = bp1l(current_node, p2, density, mx, my, hx, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, get_left (current_node), tmp);

                V_SetCmp(x, current_node + nodes_count, current_node % 2);
                tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
                V_SetCmp(b, nodes_count + current_node, tmp);
            }
        }
    }

    // Left Up Angle
    {
        Q_SetLen(a, current_node, nz_angle);
        {
            tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
            curr_nz = 0;
            Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

            tmp = ap1r(current_node, density, pulse, hx, mu, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_right (current_node), tmp);

            tmp = ap1d (mu, hy);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

            tmp = ap20 (current_node, density, pulse, hy, mx ,my);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

            tmp = ap2d(current_node, density, pulse, hy, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

            V_SetCmp(x, current_node, current_node % 2);
            tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
            V_SetCmp(b, current_node, tmp);

        }
        Q_SetLen(a, nodes_count + current_node, nz_angle);
        {
            tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
            curr_nz = 0;
            Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

            tmp = bp2r(mu, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_right (current_node), tmp);

            tmp = bp2d (current_node, p2, density, mx, my, hy, mu);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

            tmp = bp10(current_node, p2, density, mx, my, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

            tmp = bp1r(current_node, p2, density, mx, my, hx, mu);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, get_right (current_node), tmp);

            V_SetCmp(x, current_node + nodes_count, current_node % 2);
            tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
            V_SetCmp(b, nodes_count + current_node, tmp);
        }
        current_node++;
    }
    // Up side
    {
        for (; current_node < nodes_count - 1; current_node++)
        {
            Q_SetLen(a, current_node, nz_side_vert1);
            {
                tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
                curr_nz = 0;
                Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

                tmp = ap1r(current_node, density, pulse, hx, mu, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_right (current_node), tmp);

                tmp = ap1l(current_node, density, pulse, hx, mu, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

                tmp = ap1d (mu, hy);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

                tmp = ap20 (current_node, density, pulse, hy, mx ,my);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

                tmp = ap2d(current_node, density, pulse, hy, mx);
                curr_nz++;
                Q_SetEntry(a, current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                V_SetCmp(x, current_node, current_node % 2);
                tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
                V_SetCmp(b, current_node, tmp);

            }
            Q_SetLen(a, nodes_count + current_node, nz_side_vert2);
            {
                tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
                curr_nz = 0;
                Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

                tmp = bp2r(mu, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_right (current_node), tmp);

                tmp = bp2l (mu, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_left (current_node), tmp);

                tmp = bp2d (current_node, p2, density, mx, my, hy, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

                tmp = bp10(current_node, p2, density, mx, my, hx);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

                tmp = bp1r(current_node, p2, density, mx, my, hx, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, get_right (current_node), tmp);

                tmp = bp1l(current_node, p2, density, mx, my, hx, mu);
                curr_nz++;
                Q_SetEntry(a, nodes_count + current_node, curr_nz, get_left (current_node), tmp);

                V_SetCmp(x, current_node + nodes_count, current_node % 2);
                tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
                V_SetCmp(b, nodes_count + current_node, tmp);
            }
        }
    }
    // Up Right Angle
    {
        Q_SetLen(a, current_node, nz_angle);
        {
            tmp = ap10 (current_node, density, pulse, hx, hy, tau, mu, mx);
            curr_nz = 0;
            Q_SetEntry (a, current_node, curr_nz, current_node, tmp);

            tmp = ap1l(current_node, density, pulse, hx, mu, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_left (current_node), tmp);

            tmp = ap1d (mu, hy);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, get_down (current_node, mx), tmp);

            tmp = ap20 (current_node, density, pulse, hy, mx ,my);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, current_node + nodes_count, tmp);

            tmp = ap2d(current_node, density, pulse, hy, mx);
            curr_nz++;
            Q_SetEntry(a, current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

            V_SetCmp(x, current_node, current_node % 2);
            tmp = a_rhs (current_node, density, old_density, pulse, tau, GAMMA, mu, hx, hy, mx, my, nodes_count, t);
            V_SetCmp(b, current_node, tmp);

        }
        Q_SetLen(a, nodes_count + current_node, nz);
        {
            tmp = bp20(current_node, p2, density, mx, my, hx, hy, tau, mu);
            curr_nz = 0;
            Q_SetEntry(a, current_node + nodes_count, curr_nz, current_node + nodes_count, tmp);

            tmp = bp2l (mu, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_left (current_node), tmp);

            tmp = bp2d (current_node, p2, density, mx, my, hy, mu);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, nodes_count + get_down (current_node, mx), tmp);

            tmp = bp10(current_node, p2, density, mx, my, hx);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, current_node, tmp);

            tmp = bp1l(current_node, p2, density, mx, my, hx, mu);
            curr_nz++;
            Q_SetEntry(a, nodes_count + current_node, curr_nz, get_left (current_node), tmp);

            V_SetCmp(x, current_node + nodes_count, current_node % 2);
            tmp = b_rhs (current_node, pulse, old_density, density, mu, hx, hy, mx, my, nodes_count, tau, GAMMA, t);
            V_SetCmp(b, nodes_count + current_node, tmp);
        }
        current_node++;
    }
}
