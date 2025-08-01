#ifndef BOX_H
#define BOX_H

#include <math.h>

enum BoundaryType
{
    BOUNDARY_REFLECTIVE,
    BOUNDARY_TRANSMISSIVE,
    BOUNDARY_CONNECTED,
    BOUNDARY_INFLOW,
};

#define EDGE_LEFT 0
#define EDGE_RIGHT 1
#define EDGE_BOTTOM 2
#define EDGE_TOP 3

struct Box;

struct BoundaryCondition
{
    enum BoundaryType type;
    struct Box *neighbor;
    int neighbor_edge;
    double inflow_state[4]; // Only used if type == BOUNDARY_INFLOW
};

struct Box
{
    int nx, ny; // number of cells of box in x, y directions
    int size;   // total number of cells including ghost (nx + 2) x (ny + 2)

    // Structure of arrays layout: primitive variables
    double *w_1, *w_2, *w_3, *w_4; // rho, u, v, p

    // Structure of arrays layout: primitive variables
    double *u_1, *u_2, *u_3, *u_4; // rho, mom_x, mom_y, energy

    struct BoundaryCondition boundaries[4]; // L, R, B, T
};

// Allocation / deallocation
struct Box *box_create(int nx, int ny);
void box_destroy(struct Box *b);

// Index helper macro
#define INDEX(i, j, nx) ((i) + (nx + 2) * (j))

void box_set_uniform_initial_conditions(struct Box *b, double *w);

// Boundary setter
void box_set_boundary_connected(struct Box *b, int edge, struct Box *neighbor, int neighbor_edge);
void box_set_boundary_reflective(struct Box *b, int edge);
void box_set_boundary_transmissive(struct Box *b, int edge);
void box_set_boundary_inflow(struct Box *b, double *w, int edge);

// Copy connectivity for temp box
void box_copy_boundaries(struct Box *dest, const struct Box *src);

// compute ghost cells based on boundary type
void box_apply_boundary_left_edge(struct Box *b);
void box_apply_boundary_right_edge(struct Box *b);
void box_apply_boundary_bottom_edge(struct Box *b);
void box_apply_boundary_top_edge(struct Box *b);

// wrapper for x, y directions
void box_apply_boundary_x_edges(struct Box *b);
void box_apply_boundary_y_edges(struct Box *b);

// Swap function, equivalent to copy but more efficent
static inline void box_swap(struct Box **a, struct Box **b)
{
    struct Box *tmp = *a;
    *a = *b;
    *b = tmp;
}

// Transform between primitive and conserved variables
void box_prim_to_cons(struct Box *b, double gamma);
void box_cons_to_prim(struct Box *b, double gamma);

// compute the speed of sound a cell
static inline double sound_speed(double p, double rho, double gamma)
{
    return sqrt(gamma * p / rho);
}

// Compute the maximum wave speed in the box
double box_max_wave_speed(struct Box *b, double gamma);

#endif
