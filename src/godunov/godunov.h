#ifndef GODUNOV_H
#define GODUNOV_H

#include "../box/box.h"

// calculate flux
void calculate_flux_1D(double *flux, double *w_l, double *w_r, double gamma, char direction);

// x - direction sweep
void godunov_x_sweep(struct Box *curr, struct Box *temp, double dt, double dx, double gamma);

// y - direction sweep
void godunov_y_sweep(struct Box *curr, struct Box *temp, double dt, double dy, double gamma);

#endif
