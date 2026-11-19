// Mach 3 wind tunnel with a forward-facing step (Emery 1968; Woodward & Colella 1984):
// a 3 x 1 channel with a step of height 0.2 starting at x = 0.6. Three boxes; the region
// under C is solid, so there is simply no box there:
//
//   y=1 +---------+---------------------------------+
//       |    B    |               C                 |
//   0.2 +---------+---------------------------------+
//       |    A    |
//   y=0 +---------+
//      x=0       0.6                                3
//
// base = cells per 0.2 length (default 20: 300 x 100 cells)
#include "cases/cases.h"

int case_step(struct Case *c, int base, double gamma)
{
    const int b = base;
    case_begin(c, "step", 0.2 / b, 0.2 / b, gamma, 4.0, 1000000);

    struct Box *A = case_add_box(c, "A", 3 * b, b, 0, 0);          // below the inlet, left of the step
    struct Box *B = case_add_box(c, "B", 3 * b, 4 * b, 0, b);      // above A
    struct Box *C = case_add_box(c, "C", 12 * b, 4 * b, 3 * b, b); // above the step
    if (!(A && B && C))
    {
        case_destroy(c);
        return -1;
    }

    case_connect(A, EDGE_TOP, B, EDGE_BOTTOM);
    case_connect(B, EDGE_RIGHT, C, EDGE_LEFT);

    // Mach 3 inflow on the left (sound speed sqrt(1.4 * 1 / 1.4) = 1, u = 3), open outflow on the right
    double w_in[4] = {1.4, 3.0, 0.0, 1.0};
    box_set_boundary_inflow(A, w_in, EDGE_LEFT);
    box_set_boundary_inflow(B, w_in, EDGE_LEFT);
    box_set_boundary_transmissive(C, EDGE_RIGHT);

    // Walls: channel floor and ceiling, step face and step top
    box_set_boundary_reflective(A, EDGE_BOTTOM);
    box_set_boundary_reflective(A, EDGE_RIGHT);
    box_set_boundary_reflective(C, EDGE_BOTTOM);
    box_set_boundary_reflective(B, EDGE_TOP);
    box_set_boundary_reflective(C, EDGE_TOP);

    // The whole tunnel starts filled with the inflow state
    box_set_uniform_initial_conditions(A, w_in);
    box_set_uniform_initial_conditions(B, w_in);
    box_set_uniform_initial_conditions(C, w_in);

    case_finish(c);
    return 0;
}
