// Boundary conditions: the setters that describe each edge of a box, and the routines that fill
// the ghost cells before every sweep. Ghost cells hold primitive variables only; the sweeps
// reconstruct and compute fluxes from them.
#include "mesh/box.h"
#include <stdio.h>

void box_set_boundary_connected(struct Box *b, int edge, struct Box *neighbor, int neighbor_edge)
{
    if (!b || !neighbor)
    {
        fprintf(stderr, "Error: NULL Box pointer in box_set_boundary_connected.\n");
        return;
    }

    // Left must face right and bottom must face top; anything else is almost certainly a typo
    if ((edge == EDGE_LEFT && neighbor_edge != EDGE_RIGHT) ||
        (edge == EDGE_RIGHT && neighbor_edge != EDGE_LEFT) ||
        (edge == EDGE_BOTTOM && neighbor_edge != EDGE_TOP) ||
        (edge == EDGE_TOP && neighbor_edge != EDGE_BOTTOM))
    {
        fprintf(stderr, "Warning: Unusual edge pairing (%d to %d)\n", edge, neighbor_edge);
    }

    // Ghost cells are copied cell by cell, so the shared edge must have the same number of cells
    int vertical = (edge == EDGE_LEFT || edge == EDGE_RIGHT);
    if ((vertical && b->ny != neighbor->ny) || (!vertical && b->nx != neighbor->nx))
    {
        fprintf(stderr, "Error: connected edge %d has %d cells but neighbor has %d\n", edge,
                vertical ? b->ny : b->nx, vertical ? neighbor->ny : neighbor->nx);
    }

    // Only this side is set; the caller connects the neighbor back to b
    b->boundaries[edge].type = BOUNDARY_CONNECTED;
    b->boundaries[edge].neighbor = neighbor;
    b->boundaries[edge].neighbor_edge = neighbor_edge;
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
        b->boundaries[edge].inflow_state[i] = w[i]; // (rho, u, v, p)
}

void box_copy_boundaries(struct Box *dest, const struct Box *src)
{
    for (int edge = 0; edge < 4; ++edge)
    {
        dest->boundaries[edge].type = src->boundaries[edge].type;
        dest->boundaries[edge].neighbor = src->boundaries[edge].neighbor;
        dest->boundaries[edge].neighbor_edge = src->boundaries[edge].neighbor_edge;

        if (src->boundaries[edge].type == BOUNDARY_INFLOW)
        {
            for (int i = 0; i < 4; ++i)
                dest->boundaries[edge].inflow_state[i] = src->boundaries[edge].inflow_state[i];
        }
    }
}

// Cell at position k along an edge, layer g = 1..BOX_NG away from the edge, in a box of nx x ny cells.
// inside = 0 gives the ghost cell outside the edge, inside = 1 the interior cell mirrored across it.
static int edge_cell(int edge, int k, int g, int inside, int nx, int ny)
{
    int d = inside ? g : 1 - g; // offset from the first interior layer: g = 1 -> 1 (interior) or 0 (ghost)
    switch (edge)
    {
    case EDGE_LEFT:
        return INDEX(d, k, nx);
    case EDGE_RIGHT:
        return INDEX(nx + 1 - d, k, nx);
    case EDGE_BOTTOM:
        return INDEX(k, d, nx);
    default: // EDGE_TOP
        return INDEX(k, ny + 1 - d, nx);
    }
}

// Fill the BOX_NG ghost layers along one edge from the primitive variables
static void apply_edge(struct Box *b, int edge)
{
    static const char *names[4] = {"LEFT", "RIGHT", "BOTTOM", "TOP"};
    const struct BoundaryCondition *bc = &b->boundaries[edge];
    int nx = b->nx, ny = b->ny;
    int vertical = (edge == EDGE_LEFT || edge == EDGE_RIGHT);
    int len = vertical ? ny : nx;

    if (bc->type != BOUNDARY_REFLECTIVE && bc->type != BOUNDARY_TRANSMISSIVE &&
        bc->type != BOUNDARY_CONNECTED && bc->type != BOUNDARY_INFLOW)
    {
        fprintf(stderr, "Unknown boundary type on %s edge.\n", names[edge]);
        return;
    }

    // g = ghost layer (1 = next to the edge), k = position along the edge
    for (int g = 1; g <= BOX_NG; ++g)
    {
        for (int k = 1; k <= len; ++k)
        {
            int ghost = edge_cell(edge, k, g, 0, nx, ny);

            if (bc->type == BOUNDARY_REFLECTIVE)
            {
                // Mirror across the wall, flipping the normal velocity
                int mirror = edge_cell(edge, k, g, 1, nx, ny);
                b->w_1[ghost] = b->w_1[mirror];
                b->w_2[ghost] = vertical ? -b->w_2[mirror] : b->w_2[mirror];
                b->w_3[ghost] = vertical ? b->w_3[mirror] : -b->w_3[mirror];
                b->w_4[ghost] = b->w_4[mirror];
            }
            else if (bc->type == BOUNDARY_TRANSMISSIVE)
            {
                // Zero-gradient extrapolation: every ghost layer copies the edge cell
                int edge_idx = edge_cell(edge, k, 1, 1, nx, ny);
                b->w_1[ghost] = b->w_1[edge_idx];
                b->w_2[ghost] = b->w_2[edge_idx];
                b->w_3[ghost] = b->w_3[edge_idx];
                b->w_4[ghost] = b->w_4[edge_idx];
            }
            else if (bc->type == BOUNDARY_CONNECTED)
            {
                // Layer g inside the neighbor, across its opposite edge
                const struct Box *nb = bc->neighbor;
                int src = edge_cell(bc->neighbor_edge, k, g, 1, nb->nx, nb->ny);
                b->w_1[ghost] = nb->w_1[src];
                b->w_2[ghost] = nb->w_2[src];
                b->w_3[ghost] = nb->w_3[src];
                b->w_4[ghost] = nb->w_4[src];
            }
            else // BOUNDARY_INFLOW
            {
                b->w_1[ghost] = bc->inflow_state[0];
                b->w_2[ghost] = bc->inflow_state[1];
                b->w_3[ghost] = bc->inflow_state[2];
                b->w_4[ghost] = bc->inflow_state[3];
            }
        }
    }
}

void box_apply_boundary_left_edge(struct Box *b)
{
    apply_edge(b, EDGE_LEFT);
}

void box_apply_boundary_right_edge(struct Box *b)
{
    apply_edge(b, EDGE_RIGHT);
}

void box_apply_boundary_bottom_edge(struct Box *b)
{
    apply_edge(b, EDGE_BOTTOM);
}

void box_apply_boundary_top_edge(struct Box *b)
{
    apply_edge(b, EDGE_TOP);
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
