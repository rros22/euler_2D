// Sod shock tube: high pressure gas on the left half, low pressure on the right, separated by a
// diaphragm at x = 0.5 that is removed at t = 0. A single box, n cells long and 4 cells wide,
// with open ends. Useful as a quick sanity check; the exact solution is known.
#include "cases/cases.h"

int case_sod(struct Case *c, int n, int dir, double gamma)
{
    case_begin(c, "sod", 1.0 / n, 1.0 / n, gamma, 0.2, 1000000);

    int nx = dir == 0 ? n : 4, ny = dir == 0 ? 4 : n;
    struct Box *tube = case_add_box(c, "tube", nx, ny, 0, 0);
    if (!tube)
        return -1;

    for (int edge = 0; edge < 4; ++edge)
        box_set_boundary_transmissive(tube, edge);

    // Non-uniform initial data: fill the primitive variables cell by cell
    double w_left[4] = {1.0, 0.0, 0.0, 1.0};
    double w_right[4] = {0.125, 0.0, 0.0, 0.1};
    for (int j = 1; j <= ny; ++j)
        for (int i = 1; i <= nx; ++i)
        {
            int k = dir == 0 ? i : j;                      // cell number along the tube
            double x = (k - 0.5) / n;                      // its centre
            const double *w = x < 0.5 ? w_left : w_right;
            int idx = INDEX(i, j, nx);
            tube->w_1[idx] = w[0];
            tube->w_2[idx] = w[1];
            tube->w_3[idx] = w[2];
            tube->w_4[idx] = w[3];
        }

    case_finish(c);
    return 0;
}
