# euler_2D

A 2D compressible Euler solver for an ideal gas. The domain is built from rectangular **boxes** of
equal-sized cells that are glued together along their edges; every other edge is a wall, an open
boundary or an inflow.

- Finite volume, Godunov flux from the **exact Riemann solver**
- **First order** or second-order **MUSCL-Hancock** (minmod, van Leer or MC limiter)
- Dimensional splitting (X-Y-Y-X sweeps)
- Output as compact binary frames, turned into an **mp4** by a fast C renderer

## Quick start

```bash
make                                    # builds build/bin/euler_2D and build/bin/render
./build/bin/euler_2D --case=jet --scheme=muscl --t-end=1 --out-dt=0.005
./build/bin/render results/jet          # -> results/jet/rho.mp4
```

Requirements: a C99 compiler and make. `ffmpeg` (with libx264) for videos. Python with numpy and
matplotlib only if you want to read frames from Python.

## Running the solver

```bash
./build/bin/euler_2D [options]
./build/bin/euler_2D --help             # lists all options and cases
```

| Option | Meaning | Default |
|---|---|---|
| `--case=NAME` | problem to run (see below) | `jet` |
| `--base=N` | resolution, meaning depends on the case | per case |
| `--t-end=T` | end time | per case |
| `--max-steps=N` | step limit | per case |
| `--out-dt=DT` | simulation time between frames (`0` = no frames) | `t_end / 200` |
| `--out-dir=DIR` | output directory | `results/<case>` |
| `--fields=rho,u,v,p` | variables stored in each frame | `rho` |
| `--scheme=first\|muscl` | first order or MUSCL-Hancock | `first` |
| `--limiter=minmod\|vanleer\|mc` | MUSCL slope limiter | `vanleer` |
| `--cfl=C` | CFL number | `0.4` |
| `--csv` | also write per-box CSV files at the start and end | off |

### Cases

| Case | Description | `--base` | Default grid |
|---|---|---|---|
| `jet` | original setup: pressurized chamber discharging into ambient gas | base unit | 100 → 500 × 800 |
| `sod` | Sod shock tube along x | number of cells | 200 |
| `step` | Mach 3 wind tunnel with a forward-facing step | cells per 0.2 length | 20 → 300 × 100 |
| `injection` | dense subsonic jet entering lighter gas; Kelvin-Helmholtz roll-up | cells per nozzle width | 48 → 432 × 576 |

Doubling `--base` makes a run about 8 times slower (4 times the cells, half the time step).

### Output

Each run writes into its output directory:

- `frame_00000.snap`, `frame_00001.snap`, ...: one binary frame per `--out-dt` of simulation time
  (float32, interior cells only; format in `src/io/snapshot.h`). The time step is shortened when
  needed so every frame lands exactly on its time.
- `case.txt`: a readable record of the settings and box layout of the run.

The output directory is **not** cleaned. A shorter rerun leaves old higher-numbered frames behind,
which the renderer would pick up. Remove the directory first (`rm -rf results/jet`), use a new
`--out-dir`, or run `make clean-results` to delete everything in `results/`.

## Making videos

```bash
./build/bin/render results/jet                      # density, viridis -> results/jet/rho.mp4
./build/bin/render results/jet --schlieren          # shocks and vortices, gray -> schlieren_rho.mp4
./build/bin/render results/jet --field=p --log      # needs --fields=rho,p (or p) in the run
./build/bin/render results/jet --still=100 -o f.ppm # a single frame as an image, no ffmpeg needed
```

| Option | Meaning | Default |
|---|---|---|
| `--field=rho\|u\|v\|p` | variable to show | `rho` |
| `--schlieren` | numerical schlieren image of the density gradient | off |
| `--log` | colour by log10 of the field | off |
| `--range=LO,HI` | fixed colour range | min/max over all frames |
| `--cmap=viridis\|gray` | colour map | viridis (gray for schlieren) |
| `--scale=S` | output height = S × cells in y | `1` |
| `--cells` | one pixel per cell instead of the true dx/dy aspect ratio | off |
| `--stride=K` | use every K-th frame | `1` |
| `--fps=N`, `--crf=N`, `--preset=NAME` | frame rate and x264 quality (larger crf = smaller file) | 30, 26, slow |
| `-o FILE` | output file | `DIR/<field>.mp4` |

The renderer reads the frames twice: once to find a single colour range for the whole video (no
flicker), once to draw. Areas without a box (solid regions) are dark gray. Videos are small: a few
hundred frames at 432 × 576 take 0.1–0.5 MB.

### Reading frames in Python

```python
import sys; sys.path.insert(0, "tools")
import snapshot
s = snapshot.read("results/jet/frame_00100.snap")
rho = snapshot.compose(s, "rho")   # whole domain as one 2D array, NaN outside the boxes, row 0 = bottom
print(s["time"], rho.shape)
```

or from the shell: `python tools/snapshot.py results/jet/frame_00100.snap --field rho -o rho.png`.

## Writing a new case

A case is a C function that builds the boxes, connects them, sets the boundary conditions and the
initial state. Here is the complete `step` case (`src/cases/step.c`):

```c
//   y=1 +---------+---------------------------------+
//       |    B    |               C                 |
//   0.2 +---------+---------------------------------+
//       |    A    |              (solid: no box)
//   y=0 +---------+
//      x=0       0.6                                3
int case_step(struct Case *c, int base, double gamma)
{
    const int b = base;
    case_begin(c, "step", 0.2 / b, 0.2 / b, gamma, 4.0, 1000000); // name, dx, dy, gamma, t_end, max_steps

    // name, cells nx x ny, position of the lower-left cell in the global grid
    struct Box *A = case_add_box(c, "A", 3 * b, b, 0, 0);
    struct Box *B = case_add_box(c, "B", 3 * b, 4 * b, 0, b);
    struct Box *C = case_add_box(c, "C", 12 * b, 4 * b, 3 * b, b);
    if (!(A && B && C))
    {
        case_destroy(c);
        return -1;
    }

    case_connect(A, EDGE_TOP, B, EDGE_BOTTOM);
    case_connect(B, EDGE_RIGHT, C, EDGE_LEFT);

    double w_in[4] = {1.4, 3.0, 0.0, 1.0}; // rho, u, v, p
    box_set_boundary_inflow(A, w_in, EDGE_LEFT);
    box_set_boundary_inflow(B, w_in, EDGE_LEFT);
    box_set_boundary_transmissive(C, EDGE_RIGHT);
    box_set_boundary_reflective(A, EDGE_RIGHT);   // ... and the other walls

    box_set_uniform_initial_conditions(A, w_in);  // ... same for B and C
    case_finish(c);                               // computes the conserved variables
    return 0;
}
```

Then register it:

1. Declare it in `src/cases/cases.h`: `int case_step(struct Case *c, int base, double gamma);`
2. Add a line to the table in `src/cases/cases.c`:
   `{"step", case_step, 20, "Mach 3 wind tunnel with a forward-facing step"},`
   (name, function, default `--base`, description shown by `--help`).
3. `make`. New `.c` files under `src/` are picked up automatically.

### Boundary conditions

Every edge (`EDGE_LEFT`, `EDGE_RIGHT`, `EDGE_BOTTOM`, `EDGE_TOP`) has exactly one condition; edges
start as walls.

| Call | Meaning |
|---|---|
| `box_set_boundary_reflective(b, edge)` | solid wall |
| `box_set_boundary_transmissive(b, edge)` | open boundary (zero gradient) |
| `box_set_boundary_inflow(b, w, edge)` | fixed state `w = {rho, u, v, p}` outside the edge |
| `case_connect(a, edge_a, b, edge_b)` | glue two boxes together (both directions) |

### Rules

- All boxes share one cell size `dx`, `dy`.
- A connection covers a **whole edge**, and both edges must have the **same number of cells**. If a
  box borders two boxes along one side, split it. Solid regions are simply places without a box.
- `origin_i`, `origin_j` (the last two arguments of `case_add_box`) place the box in the global cell
  grid. They are only used for output and rendering, so get them right or the picture will be wrong.
- At most 8 boxes per case (`CASE_MAX_BOXES`).

### Non-uniform initial data

Fill the primitive arrays of the interior cells directly, before `case_finish`:

```c
for (int j = 1; j <= box->ny; ++j)
    for (int i = 1; i <= box->nx; ++i)
    {
        double x = (origin_i + i - 0.5) * c->dx; // cell centre
        double y = (origin_j + j - 0.5) * c->dy;
        int idx = INDEX(i, j, box->nx);
        box->w_1[idx] = ...; // rho
        box->w_2[idx] = ...; // u
        box->w_3[idx] = ...; // v
        box->w_4[idx] = ...; // p
    }
```

`src/cases/sod.c` and `src/cases/injection.c` are examples.

## Code layout

```
src/
  main.c                 command line, time loop, frame output
  mesh/box.[ch]          the Box: cell arrays, ghost layers, INDEX(), variable conversion
  mesh/boundary.c        boundary conditions (fill ghost cells)
  numerics/riemann_exact exact Riemann solver (Toro, chapter 4)
  numerics/flux          Godunov flux through one face
  numerics/reconstruction  MUSCL-Hancock face states and limiters
  numerics/sweep         x and y sweeps: faces -> fluxes -> conservative update
  solver/stepper         one X-Y-Y-X time step over all boxes
  cases/                 problem setups and the case table
  io/snapshot            binary frames (and legacy CSV)
tools/render.c           frames -> mp4 / ppm
tools/snapshot.py        frames -> numpy
tests/                   regression, physics and I/O tests
bench/                   performance benchmarks (optional)
```

## Tests

```bash
make test
```

- **Regression**: a short run of the jet geometry must be bitwise identical to the stored results in
  `tests/golden/` (first order and MUSCL). Any change to the numerics shows up here.
- **Physics**: ghost cells, exact conservation in a closed domain, x/y symmetry, convergence order on a
  smooth wave, and Toro's five shock-tube problems against the exact solution.
- **I/O and rendering**: snapshot round trip, and a small run rendered to an image and a video.

Build flags include `-ffp-contract=off` so results are reproducible across optimization levels.
