// Second-order MUSCL-Hancock reconstruction of face states, used by the sweeps when the scheme is
// RECON_MUSCL_HANCOCK. The limiter controls how sharp (and how oscillation-prone) the result is.
#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

enum Limiter
{
    LIMITER_MINMOD,   // most diffusive
    LIMITER_VAN_LEER, // smooth harmonic mean
    LIMITER_MC,       // monotonized central, least diffusive of the three
};

// Slope from the backward difference a and forward difference b (zero at extrema)
double muscl_limit(double a, double b, enum Limiter limiter);

// MUSCL-Hancock face states for n cells along a line, in primitive variables (rho, u, v, p).
// wm, w0, wp hold the left neighbour, the cell and the right neighbour as four component pointers
// indexed k = 0..n-1. normal is the velocity component along the line (1 = u, 2 = v); dt_h = dt / h.
// face_m / face_p receive the states at the cell's minus / plus face, evolved by dt / 2.
// Cells whose face states would have non-positive density or pressure fall back to zero slope.
void muscl_hancock_faces(int n, const double *const wm[4], const double *const w0[4], const double *const wp[4],
                         int normal, double dt_h, double gamma, enum Limiter limiter,
                         double *const face_m[4], double *const face_p[4]);

#endif
