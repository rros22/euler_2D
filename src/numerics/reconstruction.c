// MUSCL-Hancock reconstruction (Toro, chapter 14) in primitive variables:
//   1. a limited slope per cell from its two neighbours,
//   2. a half time step evolution of the cell centre value (the "Hancock" predictor),
//   3. face values = evolved centre value -/+ half the slope.
// The Riemann solver at each face then sees second-order accurate states in space and time.
#include "numerics/reconstruction.h"
#include <math.h>

// a = backward difference, b = forward difference. All three limiters return zero slope at
// extrema (a and b of opposite sign), which keeps the scheme free of new oscillations.
double muscl_limit(double a, double b, enum Limiter limiter)
{
    if (a * b <= 0.0)
        return 0.0;

    switch (limiter)
    {
    case LIMITER_VAN_LEER:
        return 2.0 * a * b / (a + b);
    case LIMITER_MC:
    {
        double m = fmin(fmin(2.0 * fabs(a), 2.0 * fabs(b)), 0.5 * fabs(a + b));
        return a > 0.0 ? m : -m;
    }
    default: // LIMITER_MINMOD: the smaller of the two differences
        return a > 0.0 ? fmin(a, b) : fmax(a, b);
    }
}

void muscl_hancock_faces(int n, const double *const wm[4], const double *const w0[4], const double *const wp[4],
                         int normal, double dt_h, double gamma, enum Limiter limiter,
                         double *const face_m[4], double *const face_p[4])
{
    const int tangential = 3 - normal; // the other velocity component (1 <-> 2)
    const double half = 0.5 * dt_h;    // (dt / 2) / h

    for (int k = 0; k < n; ++k)
    {
        // Cell value and limited slope (difference per cell) of each variable
        double w[4], slope[4];
        for (int c = 0; c < 4; ++c)
        {
            w[c] = w0[c][k];
            slope[c] = muscl_limit(w[c] - wm[c][k], wp[c][k] - w[c], limiter);
        }

        double rho = w[0], u_n = w[normal], p = w[3];
        double d_rho = slope[0], d_un = slope[normal], d_ut = slope[tangential], d_p = slope[3];

        // Hancock half step of the primitive equations w_t + A(w) w_x = 0, with w_x * h = slope
        double w_half[4];
        w_half[0] = rho - half * (u_n * d_rho + rho * d_un);
        w_half[normal] = u_n - half * (u_n * d_un + d_p / rho);
        w_half[tangential] = w[tangential] - half * (u_n * d_ut);
        w_half[3] = p - half * (gamma * p * d_un + u_n * d_p);

        // Values at the minus (left/bottom) and plus (right/top) faces of the cell
        double w_minus[4], w_plus[4];
        for (int c = 0; c < 4; ++c)
        {
            w_minus[c] = w_half[c] - 0.5 * slope[c];
            w_plus[c] = w_half[c] + 0.5 * slope[c];
        }

        // Positivity fallback: if either face would have rho <= 0 or p <= 0, use the first-order
        // (constant) state for this cell
        if (!(w_minus[0] > 0.0 && w_minus[3] > 0.0 && w_plus[0] > 0.0 && w_plus[3] > 0.0))
        {
            for (int c = 0; c < 4; ++c)
                w_minus[c] = w_plus[c] = w[c];
        }

        for (int c = 0; c < 4; ++c)
        {
            face_m[c][k] = w_minus[c];
            face_p[c][k] = w_plus[c];
        }
    }
}
