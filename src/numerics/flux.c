// Godunov flux: solve the 1D Riemann problem normal to the face, take the state that sits on the
// face (x/t = 0) and evaluate the Euler flux there. The tangential velocity is passively advected,
// so it is taken from the upwind side.
#include "numerics/flux.h"
#include "numerics/riemann_exact.h"
#include <stdio.h>

void calculate_flux_1D(double *flux, double *w_l, double *w_r, double gamma, char direction)
{
    double x_t = 0.0; // sample on the face itself
    double sol[3];
    const double star_p_tol = 1e-6;
    double w_l_1D[3];
    double w_r_1D[3];

    // Map (rho, u, v, p) to the 1D state (rho, normal velocity, p)
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

    riemann_solver(sol, w_l_1D, w_r_1D, x_t, gamma, star_p_tol);

    double rho = sol[0];
    double normal_vel = sol[1]; // u in x-dir, v in y-dir
    double p = sol[2];

    // Transverse velocity from the upwind state
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

    // Flux of (rho, rho u, rho v, E) through the face
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
