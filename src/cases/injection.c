// Dense subsonic jet injected through a nozzle in the top wall into a closed-top, open-bottom tank
// of lighter gas at rest. The head rolls up into a vortex pair and the beam's shear layer develops
// Kelvin-Helmholtz waves that break up into small scales (from about t = 15 at base 48; at base 32
// numerical diffusion largely suppresses them). Lengths are in nozzle widths; base = cells per
// nozzle width (default 48: 432 x 576 cells).
//
//   y=12 +--wall--+=inflow=+--wall--+
//        |        |        |        |
//        |  left  | centre | right  |   sides and bottom open
//        |        |        |        |
//   y=0  +--------+--------+--------+
//       x=0       4        5        9
#include "cases/cases.h"
#include <stdint.h>

// Deterministic pseudo-random number in [-1, 1] from an integer (reproducible symmetry breaking)
static double hash_noise(uint32_t k)
{
    k ^= k >> 16;
    k *= 0x7feb352dU;
    k ^= k >> 15;
    k *= 0x846ca68bU;
    k ^= k >> 16;
    return k / 2147483647.5 - 1.0;
}

int case_injection(struct Case *c, int base, double gamma)
{
    const int b = base;
    const int nx_side = 4 * b, ny = 12 * b;
    case_begin(c, "injection", 1.0 / b, 1.0 / b, gamma, 25.0, 1000000);

    struct Box *left = case_add_box(c, "left", nx_side, ny, 0, 0);
    struct Box *centre = case_add_box(c, "centre", b, ny, nx_side, 0);
    struct Box *right = case_add_box(c, "right", nx_side, ny, nx_side + b, 0);
    if (!(left && centre && right))
    {
        case_destroy(c);
        return -1;
    }

    case_connect(left, EDGE_RIGHT, centre, EDGE_LEFT);
    case_connect(centre, EDGE_RIGHT, right, EDGE_LEFT);

    // Pressure-matched dense subsonic jet: density ratio 2, speed 0.8 (Mach 0.96 in the jet, 0.68 in the
    // ambient gas). A low convective Mach number, 0.8 / (c_jet + c_ambient) = 0.4, keeps the shear layer unstable.
    double w_jet[4] = {2.0, 0.0, -0.8, 1.0};
    double w_ambient[4] = {1.0, 0.0, 0.0, 1.0};

    // Top: nozzle in the middle, wall elsewhere. Sides and bottom: open.
    box_set_boundary_inflow(centre, w_jet, EDGE_TOP);
    box_set_boundary_reflective(left, EDGE_TOP);
    box_set_boundary_reflective(right, EDGE_TOP);
    box_set_boundary_transmissive(left, EDGE_LEFT);
    box_set_boundary_transmissive(right, EDGE_RIGHT);
    for (int k = 0; k < c->num_boxes; ++k)
        box_set_boundary_transmissive(c->boxes[k], EDGE_BOTTOM);

    // Ambient gas at rest with 1% density noise at constant pressure, so the flow is not exactly
    // mirror symmetric. The noise depends on the global cell number, not on the box layout.
    const int nx_total = 2 * nx_side + b;
    for (int k = 0; k < c->num_boxes; ++k)
    {
        struct Box *box = c->boxes[k];
        for (int j = 1; j <= box->ny; ++j)
            for (int i = 1; i <= box->nx; ++i)
            {
                int idx = INDEX(i, j, box->nx);
                uint32_t global_cell = (uint32_t)((j - 1) * nx_total + c->origin_i[k] + i - 1);
                box->w_1[idx] = w_ambient[0] * (1.0 + 0.01 * hash_noise(global_cell));
                box->w_2[idx] = w_ambient[1];
                box->w_3[idx] = w_ambient[2];
                box->w_4[idx] = w_ambient[3];
            }
    }

    case_finish(c);
    return 0;
}
