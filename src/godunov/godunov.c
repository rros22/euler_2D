#include "godunov.h"
#include "../riemann_solver/riemann_solver.h"
#include <stdio.h>
#include <math.h>

void calculate_flux_1D(double *flux, double *w_l, double *w_r, double gamma, char direction)
{
    double x_t = 0.0;
    double sol[3];
    double star_p_tol = pow(10, -6);
    double w_l_1D[3];
    double w_r_1D[3];

    // Select the velocity component and map to 1D system for Riemann solver
    if (direction == 'x')
    {
        // x-direction: [rho, u, p]
        w_l_1D[0] = w_l[0];
        w_l_1D[1] = w_l[1]; // u
        w_l_1D[2] = w_l[3];

        w_r_1D[0] = w_r[0];
        w_r_1D[1] = w_r[1]; // u
        w_r_1D[2] = w_r[3];
    }
    else if (direction == 'y')
    {
        // y-direction: [rho, v, p]
        w_l_1D[0] = w_l[0];
        w_l_1D[1] = w_l[2]; // v
        w_l_1D[2] = w_l[3];

        w_r_1D[0] = w_r[0];
        w_r_1D[1] = w_r[2]; // v
        w_r_1D[2] = w_r[3];
    }
    else
    {
        fprintf(stderr, "Error: invalid direction '%c' in calculate_flux_1D\n", direction);
        return;
    }

    // Solve Riemann problem
    riemann_solver(sol, w_l_1D, w_r_1D, x_t, gamma, star_p_tol);

    double rho = sol[0];
    double normal_vel = sol[1]; // u in x-dir, v in y-dir
    double p = sol[2];

    // Determine transverse velocity from upwind state
    double transverse_vel;
    if (normal_vel > 0)
    {
        transverse_vel = (direction == 'x') ? w_l[2] : w_l[1]; // v or u
    }
    else
    {
        transverse_vel = (direction == 'x') ? w_r[2] : w_r[1]; // v or u
    }

    double E = 0.5 * rho * (normal_vel * normal_vel + transverse_vel * transverse_vel) + p / (gamma - 1.0);

    // Conservative flux
    flux[0] = rho * normal_vel;
    if (direction == 'x')
    {
        flux[1] = rho * normal_vel * normal_vel + p; // momentum x
        flux[2] = rho * normal_vel * transverse_vel; // momentum y
    }
    else
    {
        flux[1] = rho * normal_vel * transverse_vel; // momentum x
        flux[2] = rho * normal_vel * normal_vel + p; // momentum y
    }
    flux[3] = normal_vel * (E + p); // energy
}

static inline void extract_primitive(double *w, const struct Box *b, int idx)
{
    w[0] = b->w_1[idx]; // rho
    w[1] = b->w_2[idx]; // u
    w[2] = b->w_3[idx]; // v
    w[3] = b->w_4[idx]; // p
}

void godunov_x_sweep(struct Box *curr, struct Box *temp, double dt, double dx, double gamma)
{
    int nx = curr->nx;
    int ny = curr->ny;

    double w_l[4], w_r[4], flux[4];

    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);
            int il = INDEX(i - 1, j, nx);
            int ir = INDEX(i + 1, j, nx);

            extract_primitive(w_l, curr, il);
            extract_primitive(w_r, curr, idx);
            calculate_flux_1D(flux, w_l, w_r, gamma, 'x');

            double flux_left[4];
            flux_left[0] = flux[0];
            flux_left[1] = flux[1];
            flux_left[2] = flux[2];
            flux_left[3] = flux[3];

            extract_primitive(w_l, curr, idx);
            extract_primitive(w_r, curr, ir);
            calculate_flux_1D(flux, w_l, w_r, gamma, 'x');

            double flux_right[4];
            flux_right[0] = flux[0];
            flux_right[1] = flux[1];
            flux_right[2] = flux[2];
            flux_right[3] = flux[3];

            for (int k = 0; k < 4; ++k)
            {
                double *u = (k == 0) ? temp->u_1 : (k == 1) ? temp->u_2
                                               : (k == 2)   ? temp->u_3
                                                            : temp->u_4;

                double *u_old = (k == 0) ? curr->u_1 : (k == 1) ? curr->u_2
                                                   : (k == 2)   ? curr->u_3
                                                                : curr->u_4;

                u[idx] = u_old[idx] + (dt / dx) * (flux_left[k] - flux_right[k]);
            }
        }
    }
}

void godunov_y_sweep(struct Box *curr, struct Box *temp, double dt, double dy, double gamma)
{
    int nx = curr->nx;
    int ny = curr->ny;

    double w_l[4], w_r[4], flux[4];

    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);
            int jl = INDEX(i, j - 1, nx);
            int jr = INDEX(i, j + 1, nx);

            extract_primitive(w_l, curr, jl);
            extract_primitive(w_r, curr, idx);
            calculate_flux_1D(flux, w_l, w_r, gamma, 'y');

            double flux_bottom[4];
            flux_bottom[0] = flux[0];
            flux_bottom[1] = flux[1];
            flux_bottom[2] = flux[2];
            flux_bottom[3] = flux[3];

            extract_primitive(w_l, curr, idx);
            extract_primitive(w_r, curr, jr);
            calculate_flux_1D(flux, w_l, w_r, gamma, 'y');

            double flux_top[4];
            flux_top[0] = flux[0];
            flux_top[1] = flux[1];
            flux_top[2] = flux[2];
            flux_top[3] = flux[3];

            for (int k = 0; k < 4; ++k)
            {
                double *u = (k == 0) ? temp->u_1 : (k == 1) ? temp->u_2
                                               : (k == 2)   ? temp->u_3
                                                            : temp->u_4;

                double *u_old = (k == 0) ? curr->u_1 : (k == 1) ? curr->u_2
                                                   : (k == 2)   ? curr->u_3
                                                                : curr->u_4;

                u[idx] = u_old[idx] + (dt / dy) * (flux_bottom[k] - flux_top[k]);
            }
        }
    }
}
