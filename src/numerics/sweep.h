// Directional sweeps of the finite-volume scheme. A sweep advances every cell of a box by dt
// using fluxes in one direction only; the solver (stepper.c) combines x and y sweeps.
#ifndef SWEEP_H
#define SWEEP_H

#include "mesh/box.h"
#include "numerics/flux.h"
#include "numerics/reconstruction.h"

enum Reconstruction
{
    RECON_FIRST_ORDER,   // piecewise constant (classic Godunov)
    RECON_MUSCL_HANCOCK, // limited piecewise linear + Hancock half step, needs BOX_NG >= 2
};

struct SchemeConfig
{
    enum Reconstruction recon;
    enum Limiter limiter; // used by RECON_MUSCL_HANCOCK
};

// Sweeps read curr->w (including ghost cells) and curr->u, and write temp->u and temp->w (interior).
// temp may equal curr (in-place update). If vmax is non-NULL it receives the maximum wave speed
// of the updated cells. A NULL cfg selects the first-order scheme.
void godunov_x_sweep_cfg(struct Box *curr, struct Box *temp, double dt, double dx, double gamma,
                         const struct SchemeConfig *cfg, double *vmax);
void godunov_y_sweep_cfg(struct Box *curr, struct Box *temp, double dt, double dy, double gamma,
                         const struct SchemeConfig *cfg, double *vmax);

// First-order x and y sweeps (original API; same as the _cfg versions with cfg = NULL)
void godunov_x_sweep(struct Box *curr, struct Box *temp, double dt, double dx, double gamma);
void godunov_y_sweep(struct Box *curr, struct Box *temp, double dt, double dy, double gamma);

#endif
