#include "riemann_solver.h"
#include <math.h>
#include <stdio.h>

double sound_speed(double p, double rho, double gamma)
{
    return sqrt(gamma * p / rho);
}

double fk(double p, double *w, double A, double B, double gamma, double c)
{
    double sol;

    // primitives from state w
    double p_k = w[2];

    // solution for shock
    if (p > p_k)
    {
        sol = (p - p_k) * sqrt(A / (p + B));
        return sol;
    }
    // solution for rarefaction
    else
    {
        sol = (2 * c / (gamma - 1)) * (pow((p / p_k), ((gamma - 1) / (2 * gamma))) - 1);
        return sol;
    }
}

// f_left and f_right differentiated
double d_fk(double p, double *w, double A, double B, double gamma, double c)
{

    double sol;

    // primitives from state w
    double rho_k = w[0];
    double p_k = w[2];

    // solution for shock
    if (p > p_k)
    {
        sol = pow((A / (B + p)), 0.5) * (1 - (p - p_k) / 2 / (B + p));
        return sol;
    }

    // solution for rarefaction
    else
    {
        sol = (1 / rho_k / c) * pow((p / p_k), (-(gamma + 1) / 2 / gamma));
        return sol;
    }
}

// velocity difference
double del_u(double *w_l, double *w_r)
{
    double u_l = w_l[1];
    double u_r = w_r[1];
    return u_r - u_l;
}

// f total
double f(double p, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r)
{
    double sol_l = fk(p, w_l, A_l, B_l, gamma, c_l);
    double sol_r = fk(p, w_r, A_r, B_r, gamma, c_r);
    double sol_u = del_u(w_l, w_r);
    return sol_l + sol_r + sol_u;
}

// f total differentiated
double df(double p, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r)
{

    double sol_l = d_fk(p, w_l, A_l, B_l, gamma, c_l);
    double sol_r = d_fk(p, w_r, A_r, B_r, gamma, c_r);

    return sol_l + sol_r;
}

// tolerance check
double cha(double p_old, double p_new)
{
    return sqrt(pow((p_new - p_old), 2)) / 0.5 / (p_new + p_old);
}

// rho star
double rho_star(double p_star, double *w, double gamma)
{
    double sol;
    // primitives from state w
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

// solve star star_region
void star_region(double *store, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r, double star_p_tol)
{
    // primitives from state w_l
    double u_l = w_l[1];
    double p_l = w_l[2];

    // primitives from state w_r
    double u_r = w_r[1];
    double p_r = w_r[2];

    // P_star calculation
    double p_star;

    // Define p0 guesses
    // p_guess1
    double p_guess1 = 0.5 * (p_l + p_r);

    // p_guess2
    double num = c_l + c_r - 0.5 * (gamma - 1) * (u_r - u_l);
    double den = c_l / pow(p_l, (gamma - 1) / (2 * gamma)) + c_r / pow(p_r, (gamma - 1) / (2 * gamma));
    double p_guess2 = pow(num / den, (2 * gamma) / (gamma - 1));

    // Try guess 1 first
    double p_old = p_guess1;
    double p_new = p_old;

    while (1)
    {
        double sol_f = f(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
        double sol_df = df(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
        p_new = p_old - sol_f / sol_df;

        if (cha(p_old, p_new) < star_p_tol)
        {
            p_star = p_new;
            break;
        }

        p_old = p_new;

        if (p_new < 0)
        {
            printf("Negative star pressure computed, trying guess 2:\n");
            p_star = p_new;
            break;
        }
    }

    // Try guess 2 if it computes negative pressure
    if (p_star < 0)
    {
        p_old = p_guess2;

        while (1)
        {
            double sol_f = f(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
            double sol_df = df(p_old, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r);
            p_new = p_old - sol_f / sol_df;

            if (cha(p_old, p_new) < 0.000001)
            {
                p_star = p_new;
                break;
            }

            if (p_new < 0)
            {
                printf("Negative star pressure computed, fail!!!\n");
                p_star = p_new;
                break;
            }

            p_old = p_new;
        }
    }

    // U_star calculation
    double u_star = 0.5 * (u_l + u_r) + 0.5 * (fk(p_star, w_r, A_r, B_r, gamma, c_r) -
                                               fk(p_star, w_l, A_l, B_l, gamma, c_l));

    // rho star calculation
    double rho_star_left = rho_star(p_star, w_l, gamma);
    double rho_star_right = rho_star(p_star, w_r, gamma);

    store[0] = p_star;
    store[1] = u_star;
    store[2] = rho_star_left;
    store[3] = rho_star_right;
}

// shock speeds
double shock_speed(double *w, double p_star, double c, double gamma, bool direction)
{
    double s;
    // primitives from state w

    double u_k = w[1];
    double p_k = w[2];

    if (direction)
    {
        s = 1;
    }

    else
    {
        s = -1;
    }

    return u_k + s * c * pow((gamma + 1) / (2 * gamma) * p_star / p_k + (gamma - 1) / (2 * gamma), 0.5);
}

// rarefaction speeds
double head_speed(double *w, double c, bool direction)
{
    double s;
    double u_k = w[1];

    if (direction)
    {
        s = 1;
    }

    else
    {
        s = -1;
    }

    return u_k + s * c;
}

double tail_speed(double *w, double c, double u_star, double p_star, double gamma, bool direction)
{
    double s;

    // primitives from state w
    double p_k = w[2];

    if (direction)
    {
        s = 1;
    }

    else
    {
        s = -1;
    }

    // calculate star speed of sound
    double c_star = c * pow(p_star / p_k, (gamma - 1) / (2 * gamma));
    return u_star + s * c_star;
}

// rarefaction state

double rho_fan(double *w, double c, double gamma, double x_t, bool direction)
{
    double s;

    // primitives from state w
    double rho_k = w[0];
    double u_k = w[1];

    if (direction)
    {
        s = 1;
    }

    else
    {
        s = -1;
    }

    double rho_fan = rho_k * pow(2 / (gamma + 1) - s * (gamma - 1) / (gamma + 1) / c * (u_k - x_t), 2 / (gamma - 1));
    return rho_fan;
}

double u_fan(double *w, double c, double gamma, double x_t, bool direction)
{

    double s;

    // primitives from state w
    double u_k = w[1];

    if (direction)
    {
        s = 1;
    }

    else
    {
        s = -1;
    }

    double u_fan = 2 / (gamma + 1) * (-s * c + (gamma - 1) / 2 * u_k + x_t);
    return u_fan;
}

double p_fan(double *w, double c, double gamma, double x_t, bool direction)
{
    double s;
    // primitives from state w

    double u_k = w[1];
    double p_k = w[2];

    if (direction)
    {
        s = 1;
    }

    else
    {
        s = -1;
    }

    double p_fan = p_k * pow(2 / (gamma + 1) - s * (gamma - 1) / (gamma + 1) / c * (u_k - x_t), (2 * gamma) / (gamma - 1));

    return p_fan;
}

// sampling

// f_left

void sample_left(double *sol, double *w_l, double c_l, double gamma, double *star_region, double x_t)
{
    bool left = 0;

    // primitives
    double rho_l = w_l[0];
    double u_l = w_l[1];
    double p_l = w_l[2];

    double p_star = star_region[0];
    double u_star = star_region[1];
    double rho_star_left = star_region[2];

    double h_speed = head_speed(w_l, c_l, left);
    double t_speed = tail_speed(w_l, c_l, u_star, p_star, gamma, left);
    double s_speed = shock_speed(w_l, p_star, c_l, gamma, left);

    if (p_star > p_l)
    {
        // left shock
        if (x_t < s_speed)
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
        // left fan
        if (x_t < h_speed)
        {
            sol[0] = rho_l;
            sol[1] = u_l;
            sol[2] = p_l;
        }

        else if (x_t > t_speed)
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

// f_right

void sample_right(double *sol, double *w_r, double c_r, double gamma, double *star_region, double x_t)
{
    bool right = 1;

    // primitives
    double rho_r = w_r[0];
    double u_r = w_r[1];
    double p_r = w_r[2];

    double p_star = star_region[0];
    double u_star = star_region[1];

    double rho_star_right = star_region[3];

    double h_speed = head_speed(w_r, c_r, right);
    double t_speed = tail_speed(w_r, c_r, u_star, p_star, gamma, right);
    double s_speed = shock_speed(w_r, p_star, c_r, gamma, right);

    if (p_star > p_r)
    {
        // right shock
        if (x_t > s_speed)
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
        // right fan
        if (x_t > h_speed)
        {
            sol[0] = rho_r;
            sol[1] = u_r;
            sol[2] = p_r;
        }

        else if (x_t < t_speed)
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

// sample total

void sample(double *sol, double *w_l, double *w_r, double c_l, double c_r, double gamma, double *star_region, double x_t)
{
    double u_star = star_region[1];

    if (x_t < u_star)
    {
        sample_left(sol, w_l, c_l, gamma, star_region, x_t);
    }

    else
    {
        sample_right(sol, w_r, c_r, gamma, star_region, x_t);
    }
}

// riemann solver

void riemann_solver(double *sol, double *w_l, double *w_r, double x_t, double gamma, double star_p_tol)
{
    // calculate left and right sound speeds
    double c_l = sqrt(gamma * w_l[2] / w_l[0]);
    double c_r = sqrt(gamma * w_r[2] / w_r[0]);

    // calculate data dependant constants
    double A_l = 2 / (gamma + 1) / w_l[0];
    double A_r = 2 / (gamma + 1) / w_r[0];
    double B_l = (gamma - 1) / (gamma + 1) * w_l[2];
    double B_r = (gamma - 1) / (gamma + 1) * w_r[2];

    // calculate star star_region
    double star_region_store[4];

    star_region(star_region_store, w_l, w_r, A_r, A_l, B_r, B_l, gamma, c_l, c_r, star_p_tol);

    // sample solution
    sample(sol, w_l, w_r, c_l, c_r, gamma, star_region_store, x_t);
}
