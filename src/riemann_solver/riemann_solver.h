#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H

#include <math.h>
#include <stdbool.h>

// speed of sound
double sound_speed(double p, double rho, double gamma);

// f__left and f_right
double fk(double p, double *w, double A, double B, double gamma, double c);

// f_left and f_right differentiated
double d_fk(double p, double *w, double A, double B, double gamma, double c);

// velocity difference
double del_u(double *w_l, double *w_r);

// f total
double f(double p, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r);

// f total differentiated
double df(double p, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r);

// tolerance check
double cha(double p_old, double p_new);

// rho star
double rho_star(double p, double *w, double gamma);

// solve star region
void star_region(double *store, double *w_l, double *w_r, double A_r, double A_l, double B_r, double B_l, double gamma, double c_l, double c_r, double star_p_tol);

// shock speed
double shock_speed(double *w, double p_star, double c, double gamma, bool direction);

// rarefaction speeds
double head_speed(double *w, double c, bool direction);
double tail_speed(double *w, double c, double u_star, double p_star, double gamma, bool direction);

// rarefaction state
double rho_fan(double *w, double c, double gamma, double x_t, bool direction);
double u_fan(double *w, double c, double gamma, double x_t, bool direction);
double p_fan(double *w, double c, double gamma, double x_t, bool direction);

// sampling
// left
void sample_left(double *sol, double *w_l, double c_l, double gamma, double *star_region, double x_t);

// right
void sample_right(double *sol, double *w_r, double c_r, double gamma, double *star_region, double x_t);

// total
void sample(double *sol, double *w_l, double *w_r, double c_l, double c_r, double gamma, double *star_region, double x_t);

// riemann solver
void riemann_solver(double *sol, double *w_l, double *w_r, double x_t, double gamma, double star_p_tol);

#endif
