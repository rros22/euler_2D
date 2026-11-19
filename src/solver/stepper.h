#ifndef STEPPER_H
#define STEPPER_H

#include "mesh/box.h"
#include "numerics/sweep.h"

// Maximum wave speed over a set of boxes
double solver_max_wave_speed(struct Box **boxes, int num_boxes, double gamma);

// One X-Y-Y-X split step, updated in place; each sweep uses dt, so the solution advances by 2 * dt.
// Returns the maximum wave speed of the new state (for the next time step).
double solver_step(struct Box **boxes, int num_boxes, double dt, double dx, double dy, double gamma,
                   const struct SchemeConfig *cfg);

#endif
