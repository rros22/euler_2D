#ifndef RIEMANN_EXACT_H
#define RIEMANN_EXACT_H

// Exact solution of the 1D Riemann problem between w_l and w_r (each (rho, u, p), u normal to the
// interface), sampled at x/t = x_t: sol receives (rho, u, p). The star pressure is found by Newton
// iteration to relative tolerance star_p_tol. Vacuum is not handled.
void riemann_solver(double *sol, double *w_l, double *w_r, double x_t, double gamma, double star_p_tol);

#endif
