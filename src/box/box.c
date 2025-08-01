#include "box.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct Box *box_create(int nx, int ny)
{
    struct Box *b = (struct Box *)malloc(sizeof(struct Box));
    if (!b)
    {
        printf("Allocation failed\n");
        return NULL;
    }

    // Set the dimensions and total size of b
    b->nx = nx;
    b->ny = ny;
    b->size = (nx + 2) * (ny + 2);

    // Allocate memory for primitive variables
    b->w_1 = (double *)malloc(b->size * sizeof(double));
    b->w_2 = (double *)malloc(b->size * sizeof(double));
    b->w_3 = (double *)malloc(b->size * sizeof(double));
    b->w_4 = (double *)malloc(b->size * sizeof(double));

    // Allocate memory for conserved variables
    b->u_1 = (double *)malloc(b->size * sizeof(double));
    b->u_2 = (double *)malloc(b->size * sizeof(double));
    b->u_3 = (double *)malloc(b->size * sizeof(double));
    b->u_4 = (double *)malloc(b->size * sizeof(double));

    // If any allocation fails, destroy box
    if (!(b->w_1 && b->w_2 && b->w_3 && b->w_4 &&
          b->u_1 && b->u_2 && b->u_3 && b->u_4))
    {
        printf("Allocation failed for Box arrays\n");
        box_destroy(b);
        return NULL;
    }

    // Initialize boundaries, refelctive by default
    for (int i = 0; i < 4; ++i)
    {
        b->boundaries[i].type = BOUNDARY_REFLECTIVE;
        b->boundaries[i].neighbor = NULL;
        b->boundaries[i].neighbor_edge = -1;
    }

    return b;
}

void box_destroy(struct Box *b)
{
    if (!b)
        return;

    free(b->w_1);
    free(b->w_2);
    free(b->w_3);
    free(b->w_4);

    free(b->u_1);
    free(b->u_2);
    free(b->u_3);
    free(b->u_4);

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

void box_set_boundary_connected(struct Box *b, int edge, struct Box *neighbor, int neighbor_edge)
{
    if (!b || !neighbor)
    {
        fprintf(stderr, "Error: NULL Box pointer in box_set_boundary_connected.\n");
        return;
    }

    // Check valid edge pairings
    if ((edge == EDGE_LEFT && neighbor_edge != EDGE_RIGHT) ||
        (edge == EDGE_RIGHT && neighbor_edge != EDGE_LEFT) ||
        (edge == EDGE_BOTTOM && neighbor_edge != EDGE_TOP) ||
        (edge == EDGE_TOP && neighbor_edge != EDGE_BOTTOM))
    {
        fprintf(stderr, "Warning: Unusual edge pairing (%d to %d)\n", edge, neighbor_edge);
    }

    // Set connection on this box
    b->boundaries[edge].type = BOUNDARY_CONNECTED;
    b->boundaries[edge].neighbor = neighbor;
    b->boundaries[edge].neighbor_edge = neighbor_edge;

    // Optionally: set connection on the neighbor too
    // (Comment out this section if you're managing both sides manually)
    // neighbor->boundaries[neighbor_edge].type = BOUNDARY_CONNECTED;
    // neighbor->boundaries[neighbor_edge].neighbor = b;
    // neighbor->boundaries[neighbor_edge].neighbor_edge = edge;
}

void box_set_boundary_reflective(struct Box *b, int edge)
{
    if (!b)
        return;
    b->boundaries[edge].type = BOUNDARY_REFLECTIVE;
    b->boundaries[edge].neighbor = NULL;
    b->boundaries[edge].neighbor_edge = -1;
}

void box_set_boundary_transmissive(struct Box *b, int edge)
{
    if (!b)
        return;
    b->boundaries[edge].type = BOUNDARY_TRANSMISSIVE;
    b->boundaries[edge].neighbor = NULL;
    b->boundaries[edge].neighbor_edge = -1;
}

void box_set_boundary_inflow(struct Box *b, double *w, int edge)
{
    if (!b)
        return;
    b->boundaries[edge].type = BOUNDARY_INFLOW;
    b->boundaries[edge].neighbor = NULL;
    b->boundaries[edge].neighbor_edge = -1;

    for (int i = 0; i < 4; ++i)
        b->boundaries[edge].inflow_state[i] = w[i]; // Copy values
}

void box_copy_boundaries(struct Box *dest, const struct Box *src)
{
    for (int edge = 0; edge < 4; ++edge)
    {
        dest->boundaries[edge].type = src->boundaries[edge].type;
        dest->boundaries[edge].neighbor = src->boundaries[edge].neighbor;
        dest->boundaries[edge].neighbor_edge = src->boundaries[edge].neighbor_edge;

        // If it's an inflow, copy the inflow_state array
        if (src->boundaries[edge].type == BOUNDARY_INFLOW)
        {
            for (int i = 0; i < 4; ++i)
                dest->boundaries[edge].inflow_state[i] = src->boundaries[edge].inflow_state[i];
        }
    }
}

void box_apply_boundary_left_edge(struct Box *b)
{
    int nx = b->nx;
    int ny = b->ny;

    if (b->boundaries[EDGE_LEFT].type == BOUNDARY_REFLECTIVE)
    {
        // Reflective boundary: mirror from i = 1
        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(0, j, nx);
            int mirror = INDEX(1, j, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = -b->w_2[mirror]; // Flip u
            b->w_3[ghost] = b->w_3[mirror];  // Keep v
            b->w_4[ghost] = b->w_4[mirror];  // Keep p
        }
    }
    else if (b->boundaries[EDGE_LEFT].type == BOUNDARY_TRANSMISSIVE)
    {
        // Transmissive: copy without flipping u
        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(0, j, nx);
            int mirror = INDEX(1, j, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = b->w_2[mirror]; // No flip
            b->w_3[ghost] = b->w_3[mirror];
            b->w_4[ghost] = b->w_4[mirror];
        }
    }
    else if (b->boundaries[EDGE_LEFT].type == BOUNDARY_CONNECTED)
    {
        struct Box *neighbor = b->boundaries[EDGE_LEFT].neighbor;

        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(0, j, nx);
            int src = INDEX(neighbor->nx, j, neighbor->nx); // Right edge of neighbor

            b->w_1[ghost] = neighbor->w_1[src];
            b->w_2[ghost] = neighbor->w_2[src];
            b->w_3[ghost] = neighbor->w_3[src];
            b->w_4[ghost] = neighbor->w_4[src];
        }
    }
    else if (b->boundaries[EDGE_LEFT].type == BOUNDARY_INFLOW)
    {
        double *w = b->boundaries[EDGE_LEFT].inflow_state;
        double rho = w[0];
        double u = w[1];
        double v = w[2];
        double p = w[3];

        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(0, j, nx);
            b->w_1[ghost] = rho;
            b->w_2[ghost] = u;
            b->w_3[ghost] = v;
            b->w_4[ghost] = p;
        }
    }
    else
    {
        fprintf(stderr, "Unknown boundary type on LEFT edge.\n");
    }
}

void box_apply_boundary_right_edge(struct Box *b)
{
    int nx = b->nx;
    int ny = b->ny;

    if (b->boundaries[EDGE_RIGHT].type == BOUNDARY_REFLECTIVE)
    {
        // Reflective boundary: mirror from i = nx
        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(nx + 1, j, nx);
            int mirror = INDEX(nx, j, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = -b->w_2[mirror]; // Flip u
            b->w_3[ghost] = b->w_3[mirror];  // Keep v
            b->w_4[ghost] = b->w_4[mirror];  // Keep p
        }
    }
    else if (b->boundaries[EDGE_RIGHT].type == BOUNDARY_TRANSMISSIVE)
    {
        // Reflective boundary: mirror from i = nx
        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(nx + 1, j, nx);
            int mirror = INDEX(nx, j, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = b->w_2[mirror]; // No flip
            b->w_3[ghost] = b->w_3[mirror];
            b->w_4[ghost] = b->w_4[mirror];
        }
    }
    else if (b->boundaries[EDGE_RIGHT].type == BOUNDARY_CONNECTED)
    {
        struct Box *neighbor = b->boundaries[EDGE_RIGHT].neighbor;

        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(nx + 1, j, nx);
            int src = INDEX(1, j, neighbor->nx); // Left edge of neighbor

            b->w_1[ghost] = neighbor->w_1[src];
            b->w_2[ghost] = neighbor->w_2[src];
            b->w_3[ghost] = neighbor->w_3[src];
            b->w_4[ghost] = neighbor->w_4[src];
        }
    }
    else if (b->boundaries[EDGE_RIGHT].type == BOUNDARY_INFLOW)
    {
        double *w = b->boundaries[EDGE_RIGHT].inflow_state;
        double rho = w[0];
        double u = w[1];
        double v = w[2];
        double p = w[3];

        for (int j = 1; j <= ny; ++j)
        {
            int ghost = INDEX(nx + 1, j, nx);
            b->w_1[ghost] = rho;
            b->w_2[ghost] = u;
            b->w_3[ghost] = v;
            b->w_4[ghost] = p;
        }
    }
    else
    {
        fprintf(stderr, "Unknown boundary type on RIGHT edge.\n");
    }
}

void box_apply_boundary_bottom_edge(struct Box *b)
{
    int nx = b->nx;

    if (b->boundaries[EDGE_BOTTOM].type == BOUNDARY_REFLECTIVE)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, 0, nx);
            int mirror = INDEX(i, 1, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = b->w_2[mirror];  // keep u
            b->w_3[ghost] = -b->w_3[mirror]; // flip v
            b->w_4[ghost] = b->w_4[mirror];
        }
    }
    else if (b->boundaries[EDGE_BOTTOM].type == BOUNDARY_TRANSMISSIVE)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, 0, nx);
            int mirror = INDEX(i, 1, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = b->w_2[mirror];
            b->w_3[ghost] = b->w_3[mirror]; // No flip
            b->w_4[ghost] = b->w_4[mirror];
        }
    }
    else if (b->boundaries[EDGE_BOTTOM].type == BOUNDARY_CONNECTED)
    {
        struct Box *neighbor = b->boundaries[EDGE_BOTTOM].neighbor;

        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, 0, nx);
            int src = INDEX(i, neighbor->ny, neighbor->nx); // top edge of neighbor

            b->w_1[ghost] = neighbor->w_1[src];
            b->w_2[ghost] = neighbor->w_2[src];
            b->w_3[ghost] = neighbor->w_3[src];
            b->w_4[ghost] = neighbor->w_4[src];
        }
    }
    else if (b->boundaries[EDGE_BOTTOM].type == BOUNDARY_INFLOW)
    {
        double *w = b->boundaries[EDGE_BOTTOM].inflow_state;
        double rho = w[0];
        double u = w[1];
        double v = w[2];
        double p = w[3];

        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, 0, nx);
            b->w_1[ghost] = rho;
            b->w_2[ghost] = u;
            b->w_3[ghost] = v;
            b->w_4[ghost] = p;
        }
    }
    else
    {
        fprintf(stderr, "Unknown boundary type on BOTTOM edge.\n");
    }
}

void box_apply_boundary_top_edge(struct Box *b)
{
    int nx = b->nx;
    int ny = b->ny;

    if (b->boundaries[EDGE_TOP].type == BOUNDARY_REFLECTIVE)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, ny + 1, nx);
            int mirror = INDEX(i, ny, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = b->w_2[mirror];
            b->w_3[ghost] = -b->w_3[mirror]; // flip v
            b->w_4[ghost] = b->w_4[mirror];
        }
    }
    else if (b->boundaries[EDGE_TOP].type == BOUNDARY_TRANSMISSIVE)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, ny + 1, nx);
            int mirror = INDEX(i, ny, nx);

            b->w_1[ghost] = b->w_1[mirror];
            b->w_2[ghost] = b->w_2[mirror];
            b->w_3[ghost] = b->w_3[mirror]; // No flip
            b->w_4[ghost] = b->w_4[mirror];
        }
    }
    else if (b->boundaries[EDGE_TOP].type == BOUNDARY_CONNECTED)
    {
        struct Box *neighbor = b->boundaries[EDGE_TOP].neighbor;

        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, ny + 1, nx);
            int src = INDEX(i, 1, neighbor->nx); // bottom edge of neighbor

            b->w_1[ghost] = neighbor->w_1[src];
            b->w_2[ghost] = neighbor->w_2[src];
            b->w_3[ghost] = neighbor->w_3[src];
            b->w_4[ghost] = neighbor->w_4[src];
        }
    }
    else if (b->boundaries[EDGE_TOP].type == BOUNDARY_INFLOW)
    {
        double *w = b->boundaries[EDGE_TOP].inflow_state;
        double rho = w[0];
        double u = w[1];
        double v = w[2];
        double p = w[3];

        for (int i = 1; i <= nx; ++i)
        {
            int ghost = INDEX(i, ny + 1, nx);
            b->w_1[ghost] = rho;
            b->w_2[ghost] = u;
            b->w_3[ghost] = v;
            b->w_4[ghost] = p;
        }
    }
    else
    {
        fprintf(stderr, "Unknown boundary type on TOP edge.\n");
    }
}

void box_apply_boundary_x_edges(struct Box *b)
{
    if (!b)
        return;
    box_apply_boundary_left_edge(b);
    box_apply_boundary_right_edge(b);
}

void box_apply_boundary_y_edges(struct Box *b)
{
    if (!b)
        return;
    box_apply_boundary_bottom_edge(b);
    box_apply_boundary_top_edge(b);
}

void box_prim_to_cons(struct Box *b, double gamma)
{
    if (!b)
        return;

    int nx = b->nx;
    int ny = b->ny;

    // loop through all except for ghost cells
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
            b->u_4[idx] = 0.5 * rho * (u * u + v * v) + p / (gamma - 1.0);
        }
    }
}

void box_cons_to_prim(struct Box *b, double gamma)
{
    if (!b)
        return;

    int nx = b->nx;
    int ny = b->ny;

    // loop through all except for ghost cells
    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);

            double rho = b->u_1[idx];
            double u2 = b->u_2[idx];
            double v2 = b->u_3[idx];
            double E = b->u_4[idx];

            if (rho <= 0.0)
            {
                fprintf(stderr, "Warning: Non-positive density at (%d, %d)\n", i, j);
                rho = 1e-8; // fallback safety
            }

            double u = u2 / rho;
            double v = v2 / rho;
            double kinetic = 0.5 * rho * (u * u + v * v);
            double p = (gamma - 1.0) * (E - kinetic);

            b->w_1[idx] = rho;
            b->w_2[idx] = u;
            b->w_3[idx] = v;
            b->w_4[idx] = p;
        }
    }
}

double box_max_wave_speed(struct Box *b, double gamma)
{
    if (!b)
        return 0.0;

    int nx = b->nx;
    int ny = b->ny;

    double v_max = 0.0;

    // loop through all except for ghost cells
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
                continue;
            } // skip invalid cells

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
