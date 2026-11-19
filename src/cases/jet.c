// Original 2D jet: a high-pressure chamber on top discharges downward through an outlet column
// into ambient gas. Four boxes (cell counts for base = 100):
//
//            +---------+
//            |   top   |  chamber 100 x 200, walls on three sides
//   +--------+---------+--------+
//   |  left  | bottom  | right  |  ambient 200 x 600 | 100 x 600 | 200 x 600
//   |        |         |        |  outer edges open
//   +--------+---------+--------+
//
// The domain is the unit square split into 500 x 800 cells, so cells are 1.6 times taller than wide.
#include "cases/cases.h"

int case_jet(struct Case *c, int base, double gamma)
{
    const int nx_side = 2 * base;   // left and right boxes
    const int nx_bottom = 1 * base; // outlet column
    const int ny_bottom = 6 * base;
    const int nx_top = nx_bottom; // chamber; must match the outlet column it connects to
    const int ny_top = 2 * base;

    const double Lx = 1.0, Ly = 1.0;
    const double dx = Lx / (nx_side + nx_bottom + nx_side);
    const double dy = Ly / (ny_top + ny_bottom);
    case_begin(c, "jet", dx, dy, gamma, 10.0, 20000);

    struct Box *left = case_add_box(c, "left", nx_side, ny_bottom, 0, 0);
    struct Box *bottom = case_add_box(c, "bottom", nx_bottom, ny_bottom, nx_side, 0);
    struct Box *right = case_add_box(c, "right", nx_side, ny_bottom, nx_side + nx_bottom, 0);
    struct Box *top = case_add_box(c, "top", nx_top, ny_top, nx_side + (nx_bottom - nx_top) / 2, ny_bottom);
    if (!(left && bottom && right && top))
    {
        case_destroy(c);
        return -1;
    }

    case_connect(left, EDGE_RIGHT, bottom, EDGE_LEFT);
    case_connect(bottom, EDGE_RIGHT, right, EDGE_LEFT);
    case_connect(bottom, EDGE_TOP, top, EDGE_BOTTOM);

    // Ambient region: open to the outside on every outer edge
    box_set_boundary_transmissive(left, EDGE_LEFT);
    box_set_boundary_transmissive(left, EDGE_BOTTOM);
    box_set_boundary_transmissive(left, EDGE_TOP);
    box_set_boundary_transmissive(right, EDGE_RIGHT);
    box_set_boundary_transmissive(right, EDGE_BOTTOM);
    box_set_boundary_transmissive(right, EDGE_TOP);
    box_set_boundary_transmissive(bottom, EDGE_BOTTOM);

    // Chamber walls
    box_set_boundary_reflective(top, EDGE_TOP);
    box_set_boundary_reflective(top, EDGE_LEFT);
    box_set_boundary_reflective(top, EDGE_RIGHT);

    // Initial state (rho, u, v, p): gas at rest, chamber at 8x the density and 10x the pressure
    double w_ambient[4] = {0.125, 0.0, 0.0, 0.1};
    double w_chamber[4] = {1.0, 0.0, 0.0, 1.0};
    box_set_uniform_initial_conditions(left, w_ambient);
    box_set_uniform_initial_conditions(bottom, w_ambient);
    box_set_uniform_initial_conditions(right, w_ambient);
    box_set_uniform_initial_conditions(top, w_chamber);

    case_finish(c);
    return 0;
}
