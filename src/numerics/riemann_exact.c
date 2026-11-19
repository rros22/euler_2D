// Exact Riemann solver for the 1D Euler equations of an ideal gas (Toro, "Riemann Solvers and
// Numerical Methods for Fluid Dynamics", chapter 4). Function and variable names follow Toro:
// f_K is the pressure function of side K, A_K and B_K its data-dependent constants.
//
// 1D states are w = (rho, u, p) with u the velocity normal to the interface.
#include "numerics/riemann_exact.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

// Newton iteration cap for the star pressure (convergence normally takes a few iterations)
#define RIEMANN_MAX_ITER 100

// Toro's f_K(p): velocity jump across the left or right wave if the star pressure were p.
// Shock if p > p_K (Rankine-Hugoniot), rarefaction otherwise (isentropic).
static double fk(double p, double *w, double A, double B, double gamma, double c)
{
    double sol;

    double p_k = w[2];

    if (p > p_k)
    {
        sol = (p - p_k) * sqrt(A / (p + B));
        return sol;
    }
    else
    {
        sol = (2 * c / (gamma - 1)) * (pow((p / p_k), ((gamma - 1) / (2 * gamma))) - 1);
        return sol;
    }
}

// Derivative df_K/dp, for Newton's method
static double d_fk(double p, double *w, double A, double B, double gamma, double c)
{

    double sol;

    double rho_k = w[0];
    double p_k = w[2];

    if (p > p_k)
    {
        sol = pow((A / (B + p)), 0.5) * (1 - (p - p_k) / 2 / (B + p));
        return sol;
    }
    else
    {
        sol = (1 / rho_k / c) * pow((p / p_k), (-(gamma + 1) / 2 / gamma));
        return sol;
    }
}

static double del_u(double *w_l, double *w_r)
{
    double u_l = w_l[1];
    double u_r = w_r[1];
    return u_r - u_l;
}

// f(p) = f_L(p) + f_R(p) + (u_R - u_L); its root is the star pressure
static double f(double p, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r)
{
    double sol_l = fk(p, w_l, A_l, B_l, gamma, c_l);
    double sol_r = fk(p, w_r, A_r, B_r, gamma, c_r);
    double sol_u = del_u(w_l, w_r);
    return sol_l + sol_r + sol_u;
}

static double df(double p, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r)
{

    double sol_l = d_fk(p, w_l, A_l, B_l, gamma, c_l);
    double sol_r = d_fk(p, w_r, A_r, B_r, gamma, c_r);

    return sol_l + sol_r;
}

// Relative change between Newton iterates (Toro's convergence criterion)
static double relative_change(double p_old, double p_new)
{
    return fabs(p_new - p_old) / 0.5 / (p_new + p_old);
}

// Density on side K of the contact: shock (Rankine-Hugoniot) or rarefaction (isentropic)
static double rho_star(double p_star, double *w, double gamma)
{
    double sol;
    double rho_k = w[0];
    double p_k = w[2];

    if (p_star > p_k)
    {
        double num = (p_star / p_k) + ((gamma - 1) / (gamma + 1));
        double den = ((gamma - 1) / (gamma + 1)) * (p_star / p_k) + 1;
        sol = rho_k * (num / den);
        return sol;
    }
    else
    {
        sol = rho_k * pow((p_star / p_k), 1 / gamma);
        return sol;
    }
}

// Star region between the two nonlinear waves: store = (p*, u*, rho*_L, rho*_R).
// Newton's method from the mean pressure; if that goes negative, restart from the
// two-rarefaction estimate, which is exact when both waves are rarefactions.
static void star_region(double *store, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r, double star_p_tol)
{
    double u_l = w_l[1];
    double p_l = w_l[2];

    double u_r = w_r[1];
    double p_r = w_r[2];

    double p_star;

    // First attempt: start from the mean pressure
    double p_guess1 = 0.5 * (p_l + p_r);

    double p_old = p_guess1;
    double p_new = p_old;

    for (int iter = 0;; ++iter)
    {
        double sol_f = f(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
        double sol_df = df(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
        p_new = p_old - sol_f / sol_df;

        if (relative_change(p_old, p_new) < star_p_tol || iter >= RIEMANN_MAX_ITER)
        {
            p_star = p_new;
            break;
        }

        p_old = p_new;

        if (p_new < 0)
        {
            p_star = p_new;
            break;
        }
    }

    // Second attempt after a negative iterate: two-rarefaction initial guess
    if (p_star < 0)
    {
        double num = c_l + c_r - 0.5 * (gamma - 1) * (u_r - u_l);
        double den = c_l / pow(p_l, (gamma - 1) / (2 * gamma)) + c_r / pow(p_r, (gamma - 1) / (2 * gamma));
        double p_guess2 = pow(num / den, (2 * gamma) / (gamma - 1));

        p_old = p_guess2;

        for (int iter = 0;; ++iter)
        {
            double sol_f = f(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
            double sol_df = df(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
            p_new = p_old - sol_f / sol_df;

            if (relative_change(p_old, p_new) < 0.000001 || iter >= RIEMANN_MAX_ITER)
            {
                p_star = p_new;
                break;
            }

            if (p_new < 0)
            {
                fprintf(stderr, "Negative star pressure computed, fail!!!\n");
                p_star = p_new;
                break;
            }

            p_old = p_new;
        }
    }

    // Star velocity from the root, then the densities on each side of the contact
    double u_star = 0.5 * (u_l + u_r) + 0.5 * (fk(p_star, w_r, A_r, B_r, gamma, c_r) -
                                               fk(p_star, w_l, A_l, B_l, gamma, c_l));

    double rho_star_left = rho_star(p_star, w_l, gamma);
    double rho_star_right = rho_star(p_star, w_r, gamma);

    store[0] = p_star;
    store[1] = u_star;
    store[2] = rho_star_left;
    store[3] = rho_star_right;
}

// Wave speeds. direction = true for the right-moving wave, false for the left one.

static double shock_speed(double *w, double p_star, double c, double gamma, bool direction)
{
    double s = direction ? 1 : -1;
    double u_k = w[1];
    double p_k = w[2];

    return u_k + s * c * pow((gamma + 1) / (2 * gamma) * p_star / p_k + (gamma - 1) / (2 * gamma), 0.5);
}

// Leading edge of a rarefaction (moves into the undisturbed state)
static double head_speed(double *w, double c, bool direction)
{
    double s = direction ? 1 : -1;
    double u_k = w[1];

    return u_k + s * c;
}

// Trailing edge of a rarefaction (borders the star region)
static double tail_speed(double *w, double c, double u_star, double p_star, double gamma, bool direction)
{
    double s = direction ? 1 : -1;
    double p_k = w[2];

    double c_star = c * pow(p_star / p_k, (gamma - 1) / (2 * gamma));
    return u_star + s * c_star;
}

// State inside a rarefaction fan at x/t = x_t

static double rho_fan(double *w, double c, double gamma, double x_t, bool direction)
{
    double s = direction ? 1 : -1;
    double rho_k = w[0];
    double u_k = w[1];

    double rho_fan = rho_k * pow(2 / (gamma + 1) - s * (gamma - 1) / (gamma + 1) / c * (u_k - x_t), 2 / (gamma - 1));
    return rho_fan;
}

static double u_fan(double *w, double c, double gamma, double x_t, bool direction)
{
    double s = direction ? 1 : -1;
    double u_k = w[1];

    double u_fan = 2 / (gamma + 1) * (-s * c + (gamma - 1) / 2 * u_k + x_t);
    return u_fan;
}

static double p_fan(double *w, double c, double gamma, double x_t, bool direction)
{
    double s = direction ? 1 : -1;
    double u_k = w[1];
    double p_k = w[2];

    double p_fan = p_k * pow(2 / (gamma + 1) - s * (gamma - 1) / (gamma + 1) / c * (u_k - x_t), (2 * gamma) / (gamma - 1));

    return p_fan;
}

// Solution at x/t = x_t left of the contact: undisturbed left state, left star state,
// or inside the left rarefaction fan
static void sample_left(double *sol, double *w_l, double c_l, double gamma, double *star_region, double x_t)
{
    bool left = 0;

    double rho_l = w_l[0];
    double u_l = w_l[1];
    double p_l = w_l[2];

    double p_star = star_region[0];
    double u_star = star_region[1];
    double rho_star_left = star_region[2];

    if (p_star > p_l)
    {
        // Left shock
        if (x_t < shock_speed(w_l, p_star, c_l, gamma, left))
        {
            sol[0] = rho_l;
            sol[1] = u_l;
            sol[2] = p_l;
        }
        else
        {
            sol[0] = rho_star_left;
            sol[1] = u_star;
            sol[2] = p_star;
        }
    }
    else
    {
        // Left rarefaction
        if (x_t < head_speed(w_l, c_l, left))
        {
            sol[0] = rho_l;
            sol[1] = u_l;
            sol[2] = p_l;
        }
        else if (x_t > tail_speed(w_l, c_l, u_star, p_star, gamma, left))
        {
            sol[0] = rho_star_left;
            sol[1] = u_star;
            sol[2] = p_star;
        }
        else
        {
            sol[0] = rho_fan(w_l, c_l, gamma, x_t, left);
            sol[1] = u_fan(w_l, c_l, gamma, x_t, left);
            sol[2] = p_fan(w_l, c_l, gamma, x_t, left);
        }
    }
}

// Mirror image of sample_left for the right side of the contact
static void sample_right(double *sol, double *w_r, double c_r, double gamma, double *star_region, double x_t)
{
    bool right = 1;

    double rho_r = w_r[0];
    double u_r = w_r[1];
    double p_r = w_r[2];

    double p_star = star_region[0];
    double u_star = star_region[1];
    double rho_star_right = star_region[3];

    if (p_star > p_r)
    {
        // Right shock
        if (x_t > shock_speed(w_r, p_star, c_r, gamma, right))
        {
            sol[0] = rho_r;
            sol[1] = u_r;
            sol[2] = p_r;
        }
        else
        {
            sol[0] = rho_star_right;
            sol[1] = u_star;
            sol[2] = p_star;
        }
    }
    else
    {
        // Right rarefaction
        if (x_t > head_speed(w_r, c_r, right))
        {
            sol[0] = rho_r;
            sol[1] = u_r;
            sol[2] = p_r;
        }
        else if (x_t < tail_speed(w_r, c_r, u_star, p_star, gamma, right))
        {
            sol[0] = rho_star_right;
            sol[1] = u_star;
            sol[2] = p_star;
        }
        else
        {
            sol[0] = rho_fan(w_r, c_r, gamma, x_t, right);
            sol[1] = u_fan(w_r, c_r, gamma, x_t, right);
            sol[2] = p_fan(w_r, c_r, gamma, x_t, right);
        }
    }
}

// The contact moves at u*: points left of it see the left waves, points right of it the right ones
static void sample(double *sol, double *w_l, double *w_r, double c_l, double c_r, double gamma, double *star_region, double x_t)
{
    double u_star = star_region[1];

    if (x_t < u_star)
        sample_left(sol, w_l, c_l, gamma, star_region, x_t);
    else
        sample_right(sol, w_r, c_r, gamma, star_region, x_t);
}

void riemann_solver(double *sol, double *w_l, double *w_r, double x_t, double gamma, double star_p_tol)
{
    double c_l = sqrt(gamma * w_l[2] / w_l[0]);
    double c_r = sqrt(gamma * w_r[2] / w_r[0]);

    // Toro's data-dependent constants of the shock branch of f_K
    double A_l = 2 / (gamma + 1) / w_l[0];
    double A_r = 2 / (gamma + 1) / w_r[0];
    double B_l = (gamma - 1) / (gamma + 1) * w_l[2];
    double B_r = (gamma - 1) / (gamma + 1) * w_r[2];

    double star_region_store[4];
    star_region(star_region_store, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r, star_p_tol);

    sample(sol, w_l, w_r, c_l, c_r, gamma, star_region_store, x_t);
}
