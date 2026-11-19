// Shared test helpers: small problem setups, a time loop, conserved totals, and the golden-file
// dump/compare used by the bitwise regression test.
#include "test_util.h"
#include "solver/stepper.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void setup_clear(struct Setup *s)
{
    memset(s, 0, sizeof(*s));
}

int setup_four_box(struct Setup *s, int base, int closed, double gamma)
{
    setup_clear(s);

    const int nx_side = 2 * base;
    const int nx_bottom = base;
    const int ny_bottom = 6 * base;
    const int nx_top = nx_bottom;
    const int ny_top = 2 * base;

    s->dx = 1.0 / (nx_side + nx_bottom + nx_side);
    s->dy = 1.0 / (ny_top + ny_bottom);

    double w_low[4] = {0.125, 0.0, 0.0, 0.1};
    double w_high[4] = {1.0, 0.0, 0.0, 1.0};

    int nxs[4] = {nx_side, nx_bottom, nx_side, nx_top};
    int nys[4] = {ny_bottom, ny_bottom, ny_bottom, ny_top};
    s->n = 4;
    for (int i = 0; i < 4; ++i)
    {
        s->boxes[i] = box_create(nxs[i], nys[i]);
        if (!s->boxes[i])
            return -1;
    }

    struct Box *left = s->boxes[0], *bottom = s->boxes[1], *right = s->boxes[2], *top = s->boxes[3];

    box_set_boundary_connected(left, EDGE_RIGHT, bottom, EDGE_LEFT);
    box_set_boundary_connected(bottom, EDGE_LEFT, left, EDGE_RIGHT);
    box_set_boundary_connected(bottom, EDGE_RIGHT, right, EDGE_LEFT);
    box_set_boundary_connected(right, EDGE_LEFT, bottom, EDGE_RIGHT);
    box_set_boundary_connected(bottom, EDGE_TOP, top, EDGE_BOTTOM);
    box_set_boundary_connected(top, EDGE_BOTTOM, bottom, EDGE_TOP);

    // closed = 1 makes the domain a sealed container, so mass and energy must be conserved exactly
    void (*outer)(struct Box *, int) = closed ? box_set_boundary_reflective : box_set_boundary_transmissive;
    outer(left, EDGE_LEFT);
    outer(right, EDGE_RIGHT);
    outer(left, EDGE_BOTTOM);
    outer(right, EDGE_BOTTOM);
    outer(left, EDGE_TOP);
    outer(right, EDGE_TOP);
    outer(bottom, EDGE_BOTTOM);

    box_set_boundary_reflective(top, EDGE_TOP);
    box_set_boundary_reflective(top, EDGE_LEFT);
    box_set_boundary_reflective(top, EDGE_RIGHT);

    box_set_uniform_initial_conditions(left, w_low);
    box_set_uniform_initial_conditions(right, w_low);
    box_set_uniform_initial_conditions(bottom, w_low);
    box_set_uniform_initial_conditions(top, w_high);

    setup_finish(s, gamma);
    return 0;
}

int setup_single_box(struct Setup *s, int nx, int ny, double dx, double dy)
{
    setup_clear(s);
    s->n = 1;
    s->dx = dx;
    s->dy = dy;
    s->boxes[0] = box_create(nx, ny);
    return s->boxes[0] ? 0 : -1;
}

void setup_finish(struct Setup *s, double gamma)
{
    for (int i = 0; i < s->n; ++i)
        box_prim_to_cons(s->boxes[i], gamma);
}

void setup_destroy(struct Setup *s)
{
    for (int i = 0; i < s->n; ++i)
        box_destroy(s->boxes[i]);
    setup_clear(s);
}

double setup_run(struct Setup *s, double gamma, double cfl, double t_end, int max_steps,
                 const struct SchemeConfig *cfg)
{
    double t = 0.0;
    double vmax = solver_max_wave_speed(s->boxes, s->n, gamma);
    for (int step = 0; step < max_steps && t < t_end; ++step)
    {
        double dt = cfl * fmin(s->dx, s->dy) / vmax;
        if (t + 2.0 * dt > t_end)
            dt = 0.5 * (t_end - t);
        vmax = solver_step(s->boxes, s->n, dt, s->dx, s->dy, gamma, cfg);
        t += 2.0 * dt;
    }
    return t;
}

double setup_total(const struct Setup *s, int k)
{
    double sum = 0.0;
    for (int b = 0; b < s->n; ++b)
    {
        const struct Box *box = s->boxes[b];
        const double *u = (k == 1) ? box->u_1 : (k == 2) ? box->u_2
                                            : (k == 3)   ? box->u_3
                                                         : box->u_4;
        for (int j = 1; j <= box->ny; ++j)
            for (int i = 1; i <= box->nx; ++i)
                sum += u[INDEX(i, j, box->nx)];
    }
    return sum * s->dx * s->dy;
}

int setup_is_physical(const struct Setup *s)
{
    for (int b = 0; b < s->n; ++b)
    {
        const struct Box *box = s->boxes[b];
        for (int j = 1; j <= box->ny; ++j)
            for (int i = 1; i <= box->nx; ++i)
            {
                int idx = INDEX(i, j, box->nx);
                double rho = box->w_1[idx], p = box->w_4[idx];
                if (!isfinite(rho) || !isfinite(p) || !isfinite(box->w_2[idx]) || !isfinite(box->w_3[idx]) ||
                    rho <= 0.0 || p <= 0.0)
                    return 0;
            }
    }
    return 1;
}

static void cell_values(const struct Box *box, int idx, double *v)
{
    v[0] = box->w_1[idx];
    v[1] = box->w_2[idx];
    v[2] = box->w_3[idx];
    v[3] = box->w_4[idx];
    v[4] = box->u_1[idx];
    v[5] = box->u_2[idx];
    v[6] = box->u_3[idx];
    v[7] = box->u_4[idx];
}

// Golden file format: per box, int nx, ny, then for every interior cell 8 doubles (w_1..4, u_1..4)
int setup_dump(const struct Setup *s, const char *path)
{
    FILE *fp = fopen(path, "wb");
    if (!fp)
        return -1;
    for (int b = 0; b < s->n; ++b)
    {
        const struct Box *box = s->boxes[b];
        int dims[2] = {box->nx, box->ny};
        fwrite(dims, sizeof(int), 2, fp);
        for (int j = 1; j <= box->ny; ++j)
            for (int i = 1; i <= box->nx; ++i)
            {
                double v[8];
                cell_values(box, INDEX(i, j, box->nx), v);
                fwrite(v, sizeof(double), 8, fp);
            }
    }
    fclose(fp);
    return 0;
}

long setup_compare(const struct Setup *s, const char *path, double *max_abs_diff)
{
    FILE *fp = fopen(path, "rb");
    if (!fp)
        return -1;
    long ndiff = 0;
    *max_abs_diff = 0.0;
    for (int b = 0; b < s->n; ++b)
    {
        const struct Box *box = s->boxes[b];
        int dims[2];
        if (fread(dims, sizeof(int), 2, fp) != 2 || dims[0] != box->nx || dims[1] != box->ny)
        {
            fclose(fp);
            return -1;
        }
        for (int j = 1; j <= box->ny; ++j)
            for (int i = 1; i <= box->nx; ++i)
            {
                double v[8], g[8];
                cell_values(box, INDEX(i, j, box->nx), v);
                if (fread(g, sizeof(double), 8, fp) != 8)
                {
                    fclose(fp);
                    return -1;
                }
                for (int k = 0; k < 8; ++k)
                {
                    // Compare bit patterns, so even a 1-ulp change (or a NaN) counts as a difference
                    if (memcmp(&v[k], &g[k], sizeof(double)) != 0)
                    {
                        ++ndiff;
                        double d = fabs(v[k] - g[k]);
                        if (!(d <= *max_abs_diff))
                            *max_abs_diff = d;
                    }
                }
            }
    }
    fclose(fp);
    return ndiff;
}
