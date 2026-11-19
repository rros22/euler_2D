// Time stepping of a multi-box domain by dimensional (Strang-type) splitting: the 2D problem is
// advanced with 1D sweeps in the order X, Y, Y, X. Alternating the order cancels the leading
// splitting error, so the pair of steps is second-order accurate in time.
#include "solver/stepper.h"
#include <stddef.h>

double solver_max_wave_speed(struct Box **boxes, int num_boxes, double gamma)
{
    double vmax = 0.0;
    for (int i = 0; i < num_boxes; ++i)
    {
        double v = box_max_wave_speed(boxes[i], gamma);
        if (v > vmax)
            vmax = v;
    }
    return vmax;
}

// One x-sweep of every box. All ghost cells are filled before any box is updated, so connected
// boxes exchange data from the same time level. Returns the max wave speed if want_vmax is set.
static double x_pass(struct Box **boxes, int num_boxes, double dt, double dx, double gamma,
                     const struct SchemeConfig *cfg, int want_vmax)
{
    double vmax = 0.0;
    for (int i = 0; i < num_boxes; ++i)
        box_apply_boundary_x_edges(boxes[i]);
    for (int i = 0; i < num_boxes; ++i)
    {
        double v = 0.0;
        godunov_x_sweep_cfg(boxes[i], boxes[i], dt, dx, gamma, cfg, want_vmax ? &v : NULL);
        if (v > vmax)
            vmax = v;
    }
    return vmax;
}

// One y-sweep of every box (same pattern as x_pass)
static void y_pass(struct Box **boxes, int num_boxes, double dt, double dy, double gamma,
                   const struct SchemeConfig *cfg)
{
    for (int i = 0; i < num_boxes; ++i)
        box_apply_boundary_y_edges(boxes[i]);
    for (int i = 0; i < num_boxes; ++i)
        godunov_y_sweep_cfg(boxes[i], boxes[i], dt, dy, gamma, cfg, NULL);
}

double solver_step(struct Box **boxes, int num_boxes, double dt, double dx, double dy, double gamma,
                   const struct SchemeConfig *cfg)
{
    // The last sweep also measures the wave speed, which sets the next time step for free
    x_pass(boxes, num_boxes, dt, dx, gamma, cfg, 0);
    y_pass(boxes, num_boxes, dt, dy, gamma, cfg);
    y_pass(boxes, num_boxes, dt, dy, gamma, cfg);
    return x_pass(boxes, num_boxes, dt, dx, gamma, cfg, 1);
}
