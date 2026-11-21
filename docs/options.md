# Options

Both programs print a summary with `--help`. This page explains each option in detail.

## Solver: `euler_2D`

```bash
./build/bin/euler_2D [options]
```

Options take the form `--name=value`. An option you don't give keeps its default; several defaults
come from the case (see [Cases](cases.md)).

### Problem

| Option | Default | Meaning |
|---|---|---|
| `--case=NAME` | `jet` | Problem to run: `jet`, `sod`, `step` or `injection`. `--help` lists the cases with a one-line description. |
| `--base=N` | per case | Resolution. Its meaning depends on the case: the jet's base unit, the number of Sod cells, cells per 0.2 length for the step, cells per nozzle width for the injection. All boxes scale together, so the geometry stays the same. |

Doubling `--base` makes a run about 8 times slower: 4 times as many cells, and a time step half as
large.

### Time

Four different quantities control time. They are easy to confuse:

| Quantity | Set by | Meaning |
|---|---|---|
| End time | `--t-end=T` | The run stops when the simulation time reaches T. Without the option, the case default is used (jet 10, sod 0.2, step 4, injection 25). |
| Step limit | `--max-steps=N` | The run also stops after N steps, whichever comes first. Default from the case (jet 20000, others 1000000). |
| Frame interval | `--out-dt=DT` | A frame is written at t = 0, DT, 2 DT, ... Default `t_end / 200` (201 frames); `0` writes no frames. This does **not** affect the accuracy of the simulation. |
| Solver time step | computed | Chosen at every step as Δt = CFL · min(dx, dy) / v_max, where v_max is the fastest signal speed (flow speed plus sound speed) in the domain. One step is four sweeps (X, Y, Y, X), each by Δt, so the solution advances 2Δt per step. The step is shortened when needed to land exactly on the next frame time and on the end time. |

To run longer than the case default, pass a larger `--t-end`. If the run then stops early, the step
limit was reached; raise `--max-steps` too.

### Output

| Option | Default | Meaning |
|---|---|---|
| `--out-dir=DIR` | `results/<case>` | Directory for frames, `case.txt` and CSV files. Created if missing, never cleaned. |
| `--fields=LIST` | `rho` | Comma-separated variables stored in each frame, from `rho`, `u`, `v`, `p` (density, x velocity, y velocity, pressure). More fields mean larger frames: each field costs 4 bytes per cell. The renderer can only show fields that were stored. |
| `--csv` | off | Also write the full state of every box as CSV at the start and the end (`box_<name>_initial.csv`, `box_<name>_final.csv`) plus `box_connections.csv`, a table of every edge's boundary condition. Text output is slow and large; use it for small runs or debugging. |

### Numerical scheme

| Option | Default | Meaning |
|---|---|---|
| `--scheme=first\|muscl` | `first` | `first`: first-order Godunov (robust, diffusive). `muscl`: second-order MUSCL-Hancock, much sharper shocks and vortices for about 10% more cost. |
| `--limiter=minmod\|vanleer\|mc` | `vanleer` | Slope limiter for `muscl` (ignored for `first`). `minmod` is the most diffusive and most robust, `mc` the sharpest, `vanleer` in between. |
| `--cfl=C` | `0.4` | CFL number: fraction of a cell a wave may cross per sweep. Smaller is more robust and slower; stay below 0.5 for MUSCL. |

### Other

| Option | Meaning |
|---|---|
| `--help` | Print the options and the list of cases, then exit. |

## Renderer: `render`

```bash
./build/bin/render DIR [options]
```

`DIR` is a directory containing `frame_*.snap` files, normally `results/<case>`. The renderer reads
all frames twice: first to find one colour range for the whole video, then to draw each frame. See
[Output and rendering](output-and-rendering.md) for how images are composed.

### What to show

| Option | Default | Meaning |
|---|---|---|
| `--field=rho\|u\|v\|p` | `rho` | Variable to show. It must have been stored by the solver (`--fields`). |
| `--schlieren` | off | Numerical schlieren: brightness exp(-k · \|∇ρ\| / max\|∇ρ\|), with the maximum taken over all frames. Flat regions are white, shocks and vortex edges dark. Always uses the density. |
| `--schlieren-k=K` | `15` | Schlieren contrast. Larger K shows weaker gradients. |
| `--log` | off | Colour by log10 of the field. Useful when values span orders of magnitude. |
| `--range=LO,HI` | min/max over all frames | Fixed colour range. Values outside are clipped to the end colours. Use it to compare videos with the same scale. |
| `--cmap=viridis\|gray` | `viridis` (`gray` with `--schlieren`) | Colour map. |

### Size and quality

| Option | Default | Meaning |
|---|---|---|
| `--scale=S` | `1` | Output height = S × the number of cells in y. Use `0.5` for smaller files, `2` or more for small grids. Width and height are rounded to even numbers (needed by the video format). |
| `--cells` | off | Draw one pixel per cell. By default the image has the true physical aspect ratio, which differs when dx ≠ dy (the jet's cells are 1.6 times taller than wide). |
| `--stride=K` | `1` | Use only every K-th frame (shorter, faster video). |
| `--fps=N` | `30` | Frames per second of the video. |
| `--crf=N` | `26` | x264 quality. Larger values give smaller files; 18 is nearly lossless, 30 is still fine for smooth fields. |
| `--preset=NAME` | `slow` | x264 speed/size trade-off (`ultrafast` … `veryslow`). Slower presets give smaller files. |

### Files

| Option | Default | Meaning |
|---|---|---|
| `-o FILE` | `DIR/<field>.mp4` | Output file (`DIR/schlieren_rho.mp4` with `--schlieren`, `.ppm` with `--still`). |
| `--still=N` | off | Write frame number N as a single PPM image instead of a video. Does not need ffmpeg. |
| `--ffmpeg=PATH` | `ffmpeg` | ffmpeg executable, if it is not on your PATH. |
| `--help` | | Print the options and exit. |

## Examples

```bash
# Compare two limiters in separate directories
./build/bin/euler_2D --case=step --scheme=muscl --limiter=minmod --out-dir=results/step_minmod
./build/bin/euler_2D --case=step --scheme=muscl --limiter=mc     --out-dir=results/step_mc
./build/bin/render results/step_minmod --schlieren --scale=3
./build/bin/render results/step_mc     --schlieren --scale=3

# Pressure video with the same colour scale for several runs
./build/bin/euler_2D --case=jet --fields=rho,p --t-end=1 --out-dt=0.005
./build/bin/render results/jet --field=p --range=0.05,1

# A longer injection run at higher resolution, small video
./build/bin/euler_2D --case=injection --scheme=muscl --base=64 --t-end=40 --out-dt=0.1
./build/bin/render results/injection --scale=0.5
```
