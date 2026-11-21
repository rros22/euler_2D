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
./build/bin/euler_2D --help             # all options and cases
make test                               # regression, physics, I/O and docs checks
```

Requirements: a C99 compiler and make; `ffmpeg` (with libx264) for videos; Python with numpy and
matplotlib only to read frames from Python.

The output directory is **not** cleaned between runs: remove it first (`rm -rf results/jet`) or run
`make clean-results`.

## Cases

| Case | Description |
|---|---|
| `jet` | pressurized chamber discharging into ambient gas |
| `sod` | Sod shock tube along x |
| `step` | Mach 3 wind tunnel with a forward-facing step |
| `injection` | dense subsonic jet entering lighter gas; Kelvin-Helmholtz roll-up |

Geometry, boundary conditions, initial states and defaults of each case: [docs/cases.md](docs/cases.md).

## Documentation

| Page | Contents |
|---|---|
| [Overview and units](docs/index.md) | what the solver does, dimensionless units |
| [Getting started](docs/getting-started.md) | build, first run, first video |
| [Options](docs/options.md) | every option of the solver and the renderer |
| [Cases](docs/cases.md) | parameters of each built-in case |
| [Writing a new case](docs/new-case.md) | geometry, boundary conditions, initial data, registration |
| [Output and rendering](docs/output-and-rendering.md) | frame files, renderer, reading frames in Python |
| [Architecture](docs/architecture.md) | code map, data flow, where to change what |
| [Testing](docs/testing.md) | the tests and the golden files |
| [Maintaining](docs/maintaining.md) | conventions for code, comments and docs |
