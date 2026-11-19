#ifndef FLUX_H
#define FLUX_H

// Godunov flux through a face between primitive states w_l and w_r (rho, u, v, p),
// from the exact Riemann solution sampled at x/t = 0. direction is 'x' or 'y'.
void calculate_flux_1D(double *flux, double *w_l, double *w_r, double gamma, char direction);

#endif
