# Getting started

## Requirements

- A C99 compiler and `make` (macOS and Linux).
- `ffmpeg` with libx264, to make videos. Without it you can still write single images (`--still`).
- Optional: Python 3 with numpy and matplotlib, to read frames from Python.

## Build

```bash
make            # builds build/bin/euler_2D (solver) and build/bin/render (video tool)
make test       # optional: builds and runs all tests (takes a few seconds)
```

All build output goes into `build/`; `make clean` removes it.

## A first run (a few seconds)

The Sod shock tube is a quick check that everything works:

```bash
./build/bin/euler_2D --case=sod --scheme=muscl
```

The solver prints one line per frame and a summary:

```
case sod, scheme muscl/vanleer, cfl 0.4, t_end 0.2, frames every 0.001 -> results/sod
frame     0  step      0  t = 0.00000
...
done: 200 steps, t = 0.20000, 201 frames, 0.1 s
```

## A 2D simulation and its video

```bash
./build/bin/euler_2D --case=jet --scheme=muscl --t-end=1 --out-dt=0.005   # about 3 minutes
./build/bin/render results/jet                                            # results/jet/rho.mp4
./build/bin/render results/jet --schlieren                                # results/jet/schlieren_rho.mp4
```

For a quick preview, add `--base=40` to the solver command: a coarser grid that runs in seconds.

## Where things end up

```
results/<case>/
  frame_00000.snap ...   one binary frame per --out-dt of simulation time
  case.txt               the settings and box layout of the run
  rho.mp4                videos written by render
```

Output directories are **not** cleaned between runs. If a new run writes fewer frames than an old one,
the old higher-numbered frames stay and the renderer uses them too. Before rerunning a case, remove
its directory (`rm -rf results/jet`), or give the run its own directory with `--out-dir`.
`make clean-results` deletes all of `results/`.

## Next steps

- Every option in detail: [Options](options.md)
- What the built-in cases contain: [Cases](cases.md)
- Your own problem: [Writing a new case](new-case.md)
