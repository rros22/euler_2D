# Writing a new case

This page walks through setting up a problem from scratch. The example is a box-shaped channel
with a wider chamber, but the steps are the same for any geometry.

## 1. Sketch the geometry as boxes

The domain must be made of rectangular boxes of cells, all with the same cell size. Draw it first
and decide where the box boundaries go. Two rules shape the layout:

- **A connection covers a whole edge**, and the two connected edges must have **the same number of
  cells**. If one box borders two different boxes along one side, split it into two boxes.
- **Solid regions are simply places without a box.** Their outline becomes walls.

Example: a narrow channel (`inlet`) feeding a wider chamber. The chamber is wider than the channel,
so it is split into three boxes, and only the middle one connects to the channel:

```
   +------------+--------+------------+
   |  ch_left   | ch_mid | ch_right   |   chamber: 2 + 1 + 2 units wide, 3 units tall
   +------------+--------+------------+
                | inlet  |                channel: 1 unit wide, 2 units tall
                +--------+
```

## 2. Write the setup function

Create `src/cases/chamber.c`. A case function receives a `struct Case` to fill, the resolution
`base` and gamma, and returns 0 on success:

```c
// Narrow channel feeding a wider chamber. base = cells per unit length.
#include "cases/cases.h"

int case_chamber(struct Case *c, int base, double gamma)
{
    const int b = base;

    // Name, cell size dx and dy, gamma, default end time and step limit
    case_begin(c, "chamber", 1.0 / b, 1.0 / b, gamma, 5.0, 1000000);

    // Boxes: name, cells nx x ny, and the global cell of their lower-left corner
    struct Box *inlet = case_add_box(c, "inlet", 1 * b, 2 * b, 2 * b, 0);
    struct Box *ch_left = case_add_box(c, "ch_left", 2 * b, 3 * b, 0, 2 * b);
    struct Box *ch_mid = case_add_box(c, "ch_mid", 1 * b, 3 * b, 2 * b, 2 * b);
    struct Box *ch_right = case_add_box(c, "ch_right", 2 * b, 3 * b, 3 * b, 2 * b);
    if (!(inlet && ch_left && ch_mid && ch_right))
    {
        case_destroy(c);
        return -1;
    }

    // Connections (both directions are set by one call)
    case_connect(inlet, EDGE_TOP, ch_mid, EDGE_BOTTOM);
    case_connect(ch_left, EDGE_RIGHT, ch_mid, EDGE_LEFT);
    case_connect(ch_mid, EDGE_RIGHT, ch_right, EDGE_LEFT);

    // Every other edge is a wall unless set otherwise
    double w_in[4] = {1.0, 0.0, 1.5, 1.0}; // rho, u, v, p
    box_set_boundary_inflow(inlet, w_in, EDGE_BOTTOM);
    box_set_boundary_transmissive(ch_left, EDGE_TOP);
    box_set_boundary_transmissive(ch_mid, EDGE_TOP);
    box_set_boundary_transmissive(ch_right, EDGE_TOP);

    // Initial state: gas at rest everywhere
    double w_rest[4] = {1.0, 0.0, 0.0, 1.0};
    for (int k = 0; k < c->num_boxes; ++k)
        box_set_uniform_initial_conditions(c->boxes[k], w_rest);

    case_finish(c); // computes the conserved variables; always the last step
    return 0;
}
```

### The helpers

| Function | What it does |
|---|---|
| `case_begin(c, name, dx, dy, gamma, t_end, max_steps)` | Clears the case and sets its name, cell size and default run length. |
| `case_add_box(c, name, nx, ny, origin_i, origin_j)` | Creates a box of nx × ny cells. `origin_i`, `origin_j` give the global cell of its lower-left corner; they are only used for output and rendering. Returns NULL on failure. |
| `case_connect(a, edge_a, b, edge_b)` | Glues two boxes together along the given edges, in both directions. |
| `case_finish(c)` | Converts the initial state to conserved variables. Must be called last. |
| `case_destroy(c)` | Frees all boxes (use it when setup fails). |

### Boundary conditions

Edges are `EDGE_LEFT`, `EDGE_RIGHT`, `EDGE_BOTTOM` and `EDGE_TOP`. Each edge has exactly one
condition, and every edge starts as a wall.

| Call | Meaning |
|---|---|
| `box_set_boundary_reflective(b, edge)` | Solid wall (the default). |
| `box_set_boundary_transmissive(b, edge)` | Open boundary: waves leave with little reflection. |
| `box_set_boundary_inflow(b, w, edge)` | The state `w = {rho, u, v, p}` is held fixed just outside the edge. Point the velocity into the domain. |
| `case_connect(a, edge_a, b, edge_b)` | Connection to another box. |

### Non-uniform initial data

Instead of `box_set_uniform_initial_conditions`, fill the primitive arrays of the interior cells
before `case_finish`. Cell (i, j) of box k has its centre at:

```c
for (int k = 0; k < c->num_boxes; ++k)
{
    struct Box *box = c->boxes[k];
    for (int j = 1; j <= box->ny; ++j)
        for (int i = 1; i <= box->nx; ++i)
        {
            double x = (c->origin_i[k] + i - 0.5) * c->dx; // cell centre
            double y = (c->origin_j[k] + j - 0.5) * c->dy;
            int idx = INDEX(i, j, box->nx);
            box->w_1[idx] = ...; // rho
            box->w_2[idx] = ...; // u
            box->w_3[idx] = ...; // v
            box->w_4[idx] = ...; // p
        }
}
```

`src/cases/sod.c` (a step in x) and `src/cases/injection.c` (random noise) are working examples.

## 3. Register the case

1. Declare the function in `src/cases/cases.h`:
   ```c
   int case_chamber(struct Case *c, int base, double gamma);
   ```
2. Add a line to the `CASES[]` table in `src/cases/cases.c`: name, function, default `--base`, and
   the description shown by `--help`:
   ```c
   {"chamber", case_chamber, 32, "narrow channel feeding a wider chamber (base 32: 160x160)"},
   ```
3. Add a `## chamber` section to [docs/cases.md](cases.md) with the geometry, box table and defaults.
   `make test` fails until this section exists.
4. Build and run:
   ```bash
   make
   ./build/bin/euler_2D --case=chamber --base=16 --t-end=0.5   # quick check at low resolution
   ./build/bin/render results/chamber --still=100 -o check.ppm # look at one frame
   ```

New `.c` files in `src/` are compiled automatically; the Makefile doesn't need changes.

## 4. Check the setup

Look at the first frame (`--still=0`) to confirm the boxes are where you expect: boxes drawn in the
wrong place mean wrong origins. Then check the printout of the first steps and a few frames of the
flow.

## Common mistakes

| Symptom | Cause |
|---|---|
| `Error: connected edge ... has N cells but neighbor has M` | The two connected edges have different cell counts. Change the box sizes or split a box. |
| Boxes overlap or leave gaps in the image | `origin_i` / `origin_j` don't match the box sizes. |
| `Warning: Non-positive density` from the first step on | `case_finish` was not called (the conserved variables are still zero), or an initial state has zero or negative density or pressure. |
| `Warning: Unusual edge pairing` | Edges connected the wrong way round (left must face right, bottom must face top). |
| Strong reflections at an open boundary | Transmissive boundaries are not perfectly non-reflecting; move them further away. |
| `Error: case ... has more than 8 boxes` | `CASE_MAX_BOXES` in `cases.h` limits a case to 8 boxes; raise it if needed. |
