// One-directional finite-volume sweeps. For each line of cells along the sweep direction:
//   1. face states: cell averages (first order) or MUSCL-Hancock reconstruction,
//   2. flux at every face from the exact Riemann solver (flux.c),
//   3. conservative update  u_i += dt/h * (F_{i-1/2} - F_{i+1/2}),
//   4. conversion of the updated cells back to primitive variables, tracking the wave speed.
// The x-sweep works row by row; the y-sweep walks up the rows of a strip of columns.
#include "numerics/sweep.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Column strip width of the y-sweep (cells)
#ifndef GODUNOV_Y_STRIP
#define GODUNOV_Y_STRIP 256
#endif

// Line views: the four component arrays of a box, offset to the first cell of interest
static void prim_row(const struct Box *b, int idx, const double *w[4])
{
    w[0] = b->w_1 + idx;
    w[1] = b->w_2 + idx;
    w[2] = b->w_3 + idx;
    w[3] = b->w_4 + idx;
}

static void cons_row(const struct Box *b, int idx, double *u[4])
{
    u[0] = b->u_1 + idx;
    u[1] = b->u_2 + idx;
    u[2] = b->u_3 + idx;
    u[3] = b->u_4 + idx;
}

// F[.][k] = flux between states wl[.][k] and wr[.][k], k = 0..n-1
static void flux_line(int n, const double *const wl[4], const double *const wr[4], char direction,
                      double gamma, double *const F[4])
{
    for (int k = 0; k < n; ++k)
    {
        double l[4] = {wl[0][k], wl[1][k], wl[2][k], wl[3][k]};
        double r[4] = {wr[0][k], wr[1][k], wr[2][k], wr[3][k]};
        double flux[4];
        calculate_flux_1D(flux, l, r, gamma, direction);
        F[0][k] = flux[0];
        F[1][k] = flux[1];
        F[2][k] = flux[2];
        F[3][k] = flux[3];
    }
}

// u_new = u_old + dt/h * (F_minus - F_plus) for n cells; u_new may alias u_old
static void update_line(int n, double *const u_old[4], double *const u_new[4],
                        double *const F_minus[4], double *const F_plus[4], double dt_h)
{
    for (int c = 0; c < 4; ++c)
    {
        const double *uo = u_old[c], *fm = F_minus[c], *fp = F_plus[c];
        double *un = u_new[c];
        for (int k = 0; k < n; ++k)
            un[k] = uo[k] + dt_h * (fm[k] - fp[k]);
    }
}

// Cells i0..i0+n-1 of row j of b: conserved to primitive, tracking the maximum wave speed
static void convert_row(struct Box *b, int j, int i0, int n, double gamma, double *vmax)
{
    int nx = b->nx;
    double v_max = vmax ? *vmax : 0.0;
    for (int i = i0; i < i0 + n; ++i)
    {
        int idx = INDEX(i, j, nx);
        double w[4];
        if (!cons_to_prim_cell(b->u_1[idx], b->u_2[idx], b->u_3[idx], b->u_4[idx], gamma, w))
            fprintf(stderr, "Warning: Non-positive density at (%d, %d)\n", i, j);
        b->w_1[idx] = w[0];
        b->w_2[idx] = w[1];
        b->w_3[idx] = w[2];
        b->w_4[idx] = w[3];
        if (vmax)
        {
            double s = cell_max_wave_speed(w[0], w[1], w[2], w[3], gamma);
            if (v_max < s)
                v_max = s;
        }
    }
    if (vmax)
        *vmax = v_max;
}

// count work buffers of n cells x 4 components in one allocation; lines[l][c] is component c of buffer l
static double *alloc_lines(int count, int n, double *lines[][4])
{
    double *mem = (double *)malloc((size_t)count * 4 * n * sizeof(double));
    if (!mem)
        return NULL;
    for (int l = 0; l < count; ++l)
        for (int c = 0; c < 4; ++c)
            lines[l][c] = mem + (size_t)(4 * l + c) * n;
    return mem;
}

// Face states of n cells starting at idx0, whose neighbours along the sweep are idx -/+ stride.
// First order: both faces are the cell average (pointers into curr->w, nothing computed).
// MUSCL-Hancock: computed into scratch and returned as pointers into it.
static void line_faces(const struct Box *b, int idx0, int n, int stride, int normal, double dt_h, double gamma,
                       const struct SchemeConfig *cfg, double *const scratch_m[4], double *const scratch_p[4],
                       const double *face_m[4], const double *face_p[4])
{
    if (!cfg || cfg->recon == RECON_FIRST_ORDER)
    {
        prim_row(b, idx0, face_m);
        prim_row(b, idx0, face_p);
        return;
    }

    const double *wm[4], *w0[4], *wp[4];
    prim_row(b, idx0 - stride, wm);
    prim_row(b, idx0, w0);
    prim_row(b, idx0 + stride, wp);
    muscl_hancock_faces(n, wm, w0, wp, normal, dt_h, gamma, cfg->limiter, scratch_m, scratch_p);
    for (int c = 0; c < 4; ++c)
    {
        face_m[c] = scratch_m[c];
        face_p[c] = scratch_p[c];
    }
}

void godunov_x_sweep_cfg(struct Box *curr, struct Box *temp, double dt, double dx, double gamma,
                         const struct SchemeConfig *cfg, double *vmax)
{
    int nx = curr->nx;
    int ny = curr->ny;
    const double dt_dx = dt / dx;

    // flux[.][i] is the flux through the face between cells i and i+1, i = 0..nx.
    // faces_m / faces_p hold the minus/plus face states of cells 0..nx+1 of the current row.
    double *lines[3][4];
    double *mem = alloc_lines(3, nx + 2, lines);
    if (!mem)
    {
        fprintf(stderr, "Error: allocation failed in godunov_x_sweep\n");
        return;
    }
    double *const *flux = lines[0], *const *faces_m = lines[1], *const *faces_p = lines[2];
    if (vmax)
        *vmax = 0.0;

    for (int j = 1; j <= ny; ++j)
    {
        // Reads row j only, so the row can be overwritten afterwards even in place
        const double *face_m[4], *face_p[4], *wl[4], *wr[4];
        line_faces(curr, INDEX(0, j, nx), nx + 2, 1, 1, dt_dx, gamma, cfg, faces_m, faces_p, face_m, face_p);
        for (int c = 0; c < 4; ++c)
        {
            wl[c] = face_p[c];     // plus face of cell i
            wr[c] = face_m[c] + 1; // minus face of cell i+1
        }
        flux_line(nx + 1, wl, wr, 'x', gamma, flux);

        // Cell i (1..nx) gets flux[i-1] through its left face and flux[i] through its right face
        double *u_old[4], *u_new[4], *flux_right[4];
        cons_row(curr, INDEX(1, j, nx), u_old);
        cons_row(temp, INDEX(1, j, nx), u_new);
        for (int c = 0; c < 4; ++c)
            flux_right[c] = flux[c] + 1;
        update_line(nx, u_old, u_new, flux, flux_right, dt_dx);
        convert_row(temp, j, 1, nx, gamma, vmax);
    }

    free(mem);
}

void godunov_y_sweep_cfg(struct Box *curr, struct Box *temp, double dt, double dy, double gamma,
                         const struct SchemeConfig *cfg, double *vmax)
{
    int nx = curr->nx;
    int ny = curr->ny;
    const double dt_dy = dt / dy;
    const int stride = INDEX(0, 1, nx) - INDEX(0, 0, nx);
    const int strip = nx < GODUNOV_Y_STRIP ? nx : GODUNOV_Y_STRIP;

    // Columns are independent in a y-sweep, so the box is swept in strips of columns whose
    // working set (rows j..j+2 plus the line buffers below) stays in L1 even for very wide boxes.
    // Rolling face fluxes: F_minus below row j, F_plus above row j.
    // Two slots of minus/plus face states: row j ("cur") and row j+1 ("next").
    double *lines[6][4];
    double *mem = alloc_lines(6, strip, lines);
    if (!mem)
    {
        fprintf(stderr, "Error: allocation failed in godunov_y_sweep\n");
        return;
    }
    if (vmax)
        *vmax = 0.0;

    for (int i0 = 1; i0 <= nx; i0 += strip)
    {
        const int n = (nx + 1 - i0 < strip) ? nx + 1 - i0 : strip;
        double *const *F_minus = lines[0], *const *F_plus = lines[1];
        double *const *faces_m[2] = {lines[2], lines[4]}, *const *faces_p[2] = {lines[3], lines[5]};
        const double *cur_m[4], *cur_p[4], *next_m[4], *next_p[4];

        // Prime the loop: faces of ghost row 0 and row 1, and the flux through the bottom face of row 1
        int slot = 0;
        line_faces(curr, INDEX(i0, 0, nx), n, stride, 2, dt_dy, gamma, cfg, faces_m[slot], faces_p[slot], cur_m,
                   cur_p);
        slot ^= 1;
        line_faces(curr, INDEX(i0, 1, nx), n, stride, 2, dt_dy, gamma, cfg, faces_m[slot], faces_p[slot], next_m,
                   next_p);
        flux_line(n, cur_p, next_m, 'y', gamma, F_minus);

        for (int j = 1; j <= ny; ++j)
        {
            for (int c = 0; c < 4; ++c)
            {
                cur_m[c] = next_m[c];
                cur_p[c] = next_p[c];
            }
            // Face states of row j+1 read rows j..j+2; afterwards row j is no longer read, so the
            // update below may overwrite it even in place
            slot ^= 1;
            line_faces(curr, INDEX(i0, j + 1, nx), n, stride, 2, dt_dy, gamma, cfg, faces_m[slot], faces_p[slot],
                       next_m, next_p);
            flux_line(n, cur_p, next_m, 'y', gamma, F_plus);

            double *u_old[4], *u_new[4];
            cons_row(curr, INDEX(i0, j, nx), u_old);
            cons_row(temp, INDEX(i0, j, nx), u_new);
            update_line(n, u_old, u_new, F_minus, F_plus, dt_dy);
            convert_row(temp, j, i0, n, gamma, vmax);

            // The top face of row j is the bottom face of row j+1
            double *const *swap = F_minus;
            F_minus = F_plus;
            F_plus = swap;
        }
    }

    free(mem);
}

// Original first-order entry points, kept for API compatibility
void godunov_x_sweep(struct Box *curr, struct Box *temp, double dt, double dx, double gamma)
{
    godunov_x_sweep_cfg(curr, temp, dt, dx, gamma, NULL, NULL);
}

void godunov_y_sweep(struct Box *curr, struct Box *temp, double dt, double dy, double gamma)
{
    godunov_y_sweep_cfg(curr, temp, dt, dy, gamma, NULL, NULL);
}
