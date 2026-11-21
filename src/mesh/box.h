// A Box is one rectangular block of the mesh: nx x ny cells of equal size, surrounded by BOX_NG
// layers of ghost cells. The domain of a case is a set of boxes glued together through their
// edges (BOUNDARY_CONNECTED); every other edge is a wall, an open boundary or an inflow.
//
// Each cell stores the primitive variables w = (rho, u, v, p) and the conserved variables
// u = (rho, rho u, rho v, E), each in its own array (structure of arrays), indexed with INDEX(i, j, nx).
#ifndef BOX_H
#define BOX_H

#include <math.h>

enum BoundaryType
{
    BOUNDARY_REFLECTIVE,   // solid wall: mirror the state, flip the normal velocity
    BOUNDARY_TRANSMISSIVE, // open boundary: copy the edge cell outward (zero gradient)
    BOUNDARY_CONNECTED,    // shared with a neighbouring box: copy its cells
    BOUNDARY_INFLOW,       // fixed primitive state in the ghost cells
};

#define EDGE_LEFT 0
#define EDGE_RIGHT 1
#define EDGE_BOTTOM 2
#define EDGE_TOP 3

struct Box;

struct BoundaryCondition
{
    enum BoundaryType type;
    struct Box *neighbor;   // BOUNDARY_CONNECTED only: the box on the other side
    int neighbor_edge;      // and which of its edges touches this one
    double inflow_state[4]; // BOUNDARY_INFLOW only: (rho, u, v, p)
};

struct Box
{
    int nx, ny; // number of cells of box in x, y directions
    int size;   // total number of cells including ghost (nx + 2 * BOX_NG) x (ny + 2 * BOX_NG)

    // Structure of arrays layout: primitive variables
    double *w_1, *w_2, *w_3, *w_4; // rho, u, v, p

    // Structure of arrays layout: conserved variables
    double *u_1, *u_2, *u_3, *u_4; // rho, mom_x, mom_y, energy

    // Single allocation backing the eight arrays above (internal)
    double *data;

    struct BoundaryCondition boundaries[4]; // L, R, B, T
};

// Allocate a box of nx x ny interior cells: all values zero, all edges reflective walls.
// Returns NULL after printing an error if allocation fails. The caller frees it with box_destroy.
struct Box *box_create(int nx, int ny);

// Free a box and its arrays; NULL is allowed. Boxes connected to it keep a dangling neighbor pointer.
void box_destroy(struct Box *b);

// Array alignment, page rounding and inter-array skew in bytes (see box_create)
#define BOX_ALIGN 64
#define BOX_PAGE 16384
#define BOX_ARRAY_SKEW 192

// Ghost layers on each side (2 = stencil of second-order reconstruction)
#define BOX_NG 2

// Array index of cell (i, j). Interior cells are i = 1..nx, j = 1..ny; ghost cells are
// i = 1 - BOX_NG..0 and nx + 1..nx + BOX_NG (likewise for j). Rows (fixed j) are contiguous.
#define INDEX(i, j, nx) ((i) + BOX_NG - 1 + ((nx) + 2 * BOX_NG) * ((j) + BOX_NG - 1))

// Same primitive state w = (rho, u, v, p) in every interior cell
void box_set_uniform_initial_conditions(struct Box *b, double *w);

// Boundary setters (boundary.c). A connection is one-directional: call it for both boxes.
// Connected edges must have the same number of cells; ghost cells are copied one to one.
void box_set_boundary_connected(struct Box *b, int edge, struct Box *neighbor, int neighbor_edge);
void box_set_boundary_reflective(struct Box *b, int edge);
void box_set_boundary_transmissive(struct Box *b, int edge);
void box_set_boundary_inflow(struct Box *b, double *w, int edge);

// Copy the boundary settings of one box to another of the same size
void box_copy_boundaries(struct Box *dest, const struct Box *src);

// Fill the ghost cells of one edge from the primitive variables, according to its boundary type
void box_apply_boundary_left_edge(struct Box *b);
void box_apply_boundary_right_edge(struct Box *b);
void box_apply_boundary_bottom_edge(struct Box *b);
void box_apply_boundary_top_edge(struct Box *b);

// Both edges crossed by an x-sweep (left, right) or a y-sweep (bottom, top)
void box_apply_boundary_x_edges(struct Box *b);
void box_apply_boundary_y_edges(struct Box *b);

// Swap the cell data of two equally sized boxes; same effect as a copy but O(1).
// The Box structs stay in place so that boundaries[].neighbor pointers remain valid.
static inline void box_swap(struct Box **a, struct Box **b)
{
    struct Box *x = *a, *y = *b;
    double *tmp;
#define BOX_SWAP_FIELD(f) (tmp = x->f, x->f = y->f, y->f = tmp)
    BOX_SWAP_FIELD(w_1);
    BOX_SWAP_FIELD(w_2);
    BOX_SWAP_FIELD(w_3);
    BOX_SWAP_FIELD(w_4);
    BOX_SWAP_FIELD(u_1);
    BOX_SWAP_FIELD(u_2);
    BOX_SWAP_FIELD(u_3);
    BOX_SWAP_FIELD(u_4);
    BOX_SWAP_FIELD(data);
#undef BOX_SWAP_FIELD
}

// Convert between primitive and conserved variables over the interior cells
void box_prim_to_cons(struct Box *b, double gamma);
void box_cons_to_prim(struct Box *b, double gamma);

// Speed of sound of an ideal gas
static inline double sound_speed(double p, double rho, double gamma)
{
    return sqrt(gamma * p / rho);
}

// Conserved to primitive for one cell; returns 0 if the density was non-positive and had to be clamped
static inline int cons_to_prim_cell(double rho, double mom_x, double mom_y, double E, double gamma, double *w)
{
    int ok = 1;
    if (rho <= 0.0)
    {
        rho = 1e-8; // fallback safety
        ok = 0;
    }
    double u = mom_x / rho;
    double v = mom_y / rho;
    double kinetic = 0.5 * rho * (u * u + v * v);
    w[0] = rho;
    w[1] = u;
    w[2] = v;
    w[3] = (gamma - 1.0) * (E - kinetic);
    return ok;
}

// max(|u| + c, |v| + c) for one cell; NaN-producing cells are ignored by the caller's comparison
static inline double cell_max_wave_speed(double rho, double u, double v, double p, double gamma)
{
    double c = sound_speed(p, rho, gamma);
    double speed_x = fabs(u) + c;
    double speed_y = fabs(v) + c;
    return (speed_x < speed_y) ? speed_y : speed_x;
}

// Maximum of |u| + c and |v| + c over the interior (for the CFL time step)
double box_max_wave_speed(struct Box *b, double gamma);

#endif
