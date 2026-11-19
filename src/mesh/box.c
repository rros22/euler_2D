// Box memory management and whole-box operations: allocation, initial data, conversion between
// primitive and conserved variables, and the maximum wave speed used for the time step.
// Boundary conditions live in boundary.c.
#define _POSIX_C_SOURCE 200112L
#include "mesh/box.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct Box *box_create(int nx, int ny)
{
    struct Box *b = (struct Box *)malloc(sizeof(struct Box));
    if (!b)
    {
        fprintf(stderr, "Allocation failed\n");
        return NULL;
    }

    b->nx = nx;
    b->ny = ny;
    b->size = (nx + 2 * BOX_NG) * (ny + 2 * BOX_NG); // interior plus ghost layers

    // One block for all eight arrays. Each array is rounded up to a whole page and successive arrays
    // are shifted by an extra BOX_ARRAY_SKEW bytes, so the same cell of different arrays does not map
    // to the same cache set (separate equal-size mallocs do, which slows the y-sweep on wide boxes).
    size_t array_bytes = ((size_t)b->size * sizeof(double) + BOX_PAGE - 1) / BOX_PAGE * BOX_PAGE;
    size_t stride = array_bytes + BOX_ARRAY_SKEW;
    void *block = NULL;
    if (posix_memalign(&block, BOX_ALIGN, 8 * stride) != 0)
    {
        fprintf(stderr, "Allocation failed for Box arrays\n");
        free(b);
        return NULL;
    }
    // Zero everything so unused ghost corners hold defined values
    memset(block, 0, 8 * stride);
    b->data = (double *)block;
    char *base = (char *)block;
    double **arrays[8] = {&b->w_1, &b->w_2, &b->w_3, &b->w_4, &b->u_1, &b->u_2, &b->u_3, &b->u_4};
    for (int k = 0; k < 8; ++k)
        *arrays[k] = (double *)(base + k * stride);

    // All edges start as reflective walls; cases overwrite them
    for (int edge = 0; edge < 4; ++edge)
    {
        b->boundaries[edge].type = BOUNDARY_REFLECTIVE;
        b->boundaries[edge].neighbor = NULL;
        b->boundaries[edge].neighbor_edge = -1;
    }

    return b;
}

void box_destroy(struct Box *b)
{
    if (!b)
        return;

    free(b->data);
    free(b);
}

void box_set_uniform_initial_conditions(struct Box *b, double *w)
{
    if (b == NULL || w == NULL)
        return;

    int nx = b->nx;
    int ny = b->ny;

    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);
            b->w_1[idx] = w[0]; // rho
            b->w_2[idx] = w[1]; // u
            b->w_3[idx] = w[2]; // v
            b->w_4[idx] = w[3]; // p
        }
    }
}

// Interior cells only: ghost cells are filled from the primitive variables by the boundary code
void box_prim_to_cons(struct Box *b, double gamma)
{
    if (!b)
        return;

    int nx = b->nx;
    int ny = b->ny;

    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);

            double rho = b->w_1[idx];
            double u = b->w_2[idx];
            double v = b->w_3[idx];
            double p = b->w_4[idx];

            b->u_1[idx] = rho;
            b->u_2[idx] = rho * u;
            b->u_3[idx] = rho * v;
            b->u_4[idx] = 0.5 * rho * (u * u + v * v) + p / (gamma - 1.0); // kinetic + internal energy
        }
    }
}

void box_cons_to_prim(struct Box *b, double gamma)
{
    if (!b)
        return;

    int nx = b->nx;
    int ny = b->ny;

    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);
            double w[4];

            if (!cons_to_prim_cell(b->u_1[idx], b->u_2[idx], b->u_3[idx], b->u_4[idx], gamma, w))
                fprintf(stderr, "Warning: Non-positive density at (%d, %d)\n", i, j);

            b->w_1[idx] = w[0];
            b->w_2[idx] = w[1];
            b->w_3[idx] = w[2];
            b->w_4[idx] = w[3];
        }
    }
}

// Largest signal speed max(|u| + c, |v| + c) over the interior; sets the CFL time step
double box_max_wave_speed(struct Box *b, double gamma)
{
    if (!b)
        return 0.0;

    int nx = b->nx;
    int ny = b->ny;

    double v_max = 0.0;

    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);

            double rho = b->w_1[idx];
            double u = b->w_2[idx];
            double v = b->w_3[idx];
            double p = b->w_4[idx];

            if (rho <= 0.0)
            {
                fprintf(stderr, "Max wave speed calculation warning: Non-positive density at (%d, %d)\n", i, j);
                continue; // skip invalid cells
            }

            double c = sound_speed(p, rho, gamma);

            double speed_x = fabs(u) + c;
            double speed_y = fabs(v) + c;

            if (v_max < speed_x)
                v_max = speed_x;
            if (v_max < speed_y)
                v_max = speed_y;
        }
    }

    return v_max;
}
