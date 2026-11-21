# Cases

A case is a complete problem setup: boxes, boundary conditions, initial state and default run
length. Each one lives in its own file in `src/cases/` and is listed in the table in
`src/cases/cases.c`. States are given as (ρ, u, v, p) in the units of [the index page](index.md#units).

Runtimes are rough figures from the development machine (the solver uses one core); scale them to
your hardware.

| Case | Purpose | Default grid | Default `t_end` |
|---|---|---|---|
| [jet](#jet) | the original setup: pressurized chamber discharging into ambient gas | 500 × 800 | 10 |
| [sod](#sod) | 1D shock tube, quick check against a known solution | 200 × 4 | 0.2 |
| [step](#step) | Mach 3 flow over a forward-facing step (classic 2D benchmark) | 300 × 100 | 4 |
| [injection](#injection) | dense jet entering lighter gas, Kelvin-Helmholtz roll-up | 432 × 576 | 25 |

## jet

The original problem of the code (`src/cases/jet.c`). A chamber of high-pressure gas on top
discharges downward through an outlet column into low-pressure ambient gas, producing a blast
wave and a starting jet with vortex roll-up.

```
            +---------+
            |   top   |  chamber, walls on three sides
   +--------+---------+--------+
   |  left  | bottom  | right  |  ambient gas, outer edges open
   |        |         |        |
   +--------+---------+--------+
```

| Box | Cells (nx × ny) | Origin (i, j) | Left | Right | Bottom | Top |
|---|---|---|---|---|---|---|
| left | 2b × 6b | (0, 0) | open | → bottom | open | open |
| bottom | b × 6b | (2b, 0) | → left | → right | open | → top |
| right | 2b × 6b | (3b, 0) | → bottom | open | open | open |
| top | b × 2b | (2b, 6b) | wall | wall | → bottom | wall |

"→ box" means connected to that box; "open" is transmissive (zero gradient).

| Parameter | Value |
|---|---|
| `--base` (b) | 100: a 500 × 800 grid |
| Domain | unit square, so dx = 1/(5b), dy = 1/(8b): cells are 1.6 times taller than wide |
| Ambient state (left, bottom, right) | (0.125, 0, 0, 0.1) |
| Chamber state (top) | (1, 0, 0, 1) |
| Reference scales | L = domain width, chamber density and pressure |
| `t_end`, `max_steps` | 10, 20000 |
| Runtime | about 45 ms per step at base 100; `--t-end=1` takes about 3 minutes |

With the defaults the step limit is reached at t ≈ 9.3, before `t_end`. The interesting part of the
flow happens in the first t ≈ 1; use `--t-end=1`.

## sod

Sod's shock tube (`src/cases/sod.c`): gas at high pressure on the left half and low pressure on the
right, separated at x = 0.5. The exact solution is known (a rarefaction, a contact and a shock), so it
is a quick sanity check. It is a single box of n × 4 cells with all edges open; the 2D solver runs it
as a 1D problem.

| Parameter | Value |
|---|---|
| `--base` (n) | 200 cells along x |
| Domain | [0, 1] × [0, 4/n], dx = dy = 1/n |
| Left state (x < 0.5) | (1, 0, 0, 1) |
| Right state (x > 0.5) | (0.125, 0, 0, 0.1) |
| `t_end`, `max_steps` | 0.2, 1000000 |
| Runtime | well under a second |

The tests run this problem (and Toro's other shock tubes) and compare with the exact solution.

## step

Mach 3 wind tunnel with a forward-facing step (`src/cases/step.c`), a standard test since Emery
(1968) and Woodward & Colella (1984). A bow shock forms in front of the step and reflects off the
channel walls, and a Mach stem develops near the top wall.

```
  y=1 +---------+---------------------------------+
      |    B    |               C                 |
  0.2 +---------+---------------------------------+
      |    A    |   (step: solid, no box)
  y=0 +---------+
     x=0       0.6                                3
```

| Box | Cells (nx × ny) | Origin (i, j) | Left | Right | Bottom | Top |
|---|---|---|---|---|---|---|
| A | 3b × b | (0, 0) | inflow | wall (step face) | wall | → B |
| B | 3b × 4b | (0, b) | inflow | → C | → A | wall |
| C | 12b × 4b | (3b, b) | → B | open | wall (step top) | wall |

| Parameter | Value |
|---|---|
| `--base` (b) | 20 cells per 0.2 length: a 300 × 100 grid |
| Domain | 3 × 1 channel, step of height 0.2 from x = 0.6, dx = dy = 0.2/b |
| Inflow and initial state | (1.4, 3, 0, 1): sound speed 1, so Mach 3 |
| Reference scales | L = channel height, ρ_ref and p_ref such that the inflow has ρ = 1.4, p = 1 |
| `t_end`, `max_steps` | 4, 1000000 |
| Runtime (MUSCL) | about 20 s at base 20, 3 min at base 40, 25 min at base 80 |

The sharp corner of the step produces a thin spurious layer along the step top; it gets thinner
with resolution. Woodward and Colella add a local fix at the corner, which this case does not.

## injection

A dense subsonic jet enters a tank of lighter gas at rest through a nozzle in the top wall
(`src/cases/injection.c`). The head rolls up into a vortex pair, and from about t ≈ 15 the edges of
the beam develop Kelvin-Helmholtz waves that break up into small-scale mixing.

```
   y=12 +--wall--+=inflow=+--wall--+
        |        |        |        |
        |  left  | centre | right  |   sides and bottom open
        |        |        |        |
   y=0  +--------+--------+--------+
       x=0       4        5        9
```

| Box | Cells (nx × ny) | Origin (i, j) | Left | Right | Bottom | Top |
|---|---|---|---|---|---|---|
| left | 4b × 12b | (0, 0) | open | → centre | open | wall |
| centre | b × 12b | (4b, 0) | → left | → right | open | inflow (jet) |
| right | 4b × 12b | (5b, 0) | → centre | open | open | wall |

| Parameter | Value |
|---|---|
| `--base` (b) | 48 cells per nozzle width: a 432 × 576 grid |
| Domain | 9 × 12 nozzle widths, dx = dy = 1/b |
| Jet (inflow) state | (2, 0, -0.8, 1): twice the tank density, moving down at 0.8, same pressure |
| Tank state | (1, 0, 0, 1), with 1% random density noise (fixed seed) to break the left-right symmetry |
| Reference scales | L = nozzle width, tank density and pressure |
| `t_end`, `max_steps` | 25, 1000000 |
| Runtime (MUSCL) | about 1.5 min at base 32, 5 min at base 48, 12 min at base 64, 40 min at base 96 |

The jet speed is low compared with the sound speeds (convective Mach number 0.4); faster jets keep
a laminar beam in 2D. At base 32 numerical diffusion suppresses most of the roll-up; use 48 or more.
The jet head reaches the bottom of the tank at t ≈ 30 and then leaves through the open boundary.
The inflow fixes all four variables, which for a subsonic jet partly reflects sound waves at the
nozzle; it doesn't visibly affect the flow.
