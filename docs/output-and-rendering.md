# Output and rendering

The solver writes compact binary frames during the run; videos and images are made afterwards
from those frames. Rendering never needs to rerun the simulation, so you can try different fields,
colour maps and sizes on the same data.

```
euler_2D  ->  results/<case>/frame_*.snap  ->  render     ->  .mp4 / .ppm
                                           ->  snapshot.py ->  numpy arrays
```

## Files written by the solver

| File | Contents |
|---|---|
| `frame_NNNNN.snap` | One frame: the stored fields of every interior cell at one time. Written at t = 0, `--out-dt`, 2 `--out-dt`, ... and at the final time. |
| `case.txt` | Plain-text record of the run settings (scheme, CFL, dx, dy, end time, fields) and the box layout. Written once at the start; no tool reads it. |
| `box_<name>_initial.csv`, `box_<name>_final.csv`, `box_connections.csv` | Only with `--csv`. |

A frame stores 4 bytes per cell per field: at the default jet resolution (400,000 cells) a
density-only frame is 1.6 MB, so 200 frames take about 320 MB. Store only the fields you need.

### Frame format

All numbers are little-endian. Frames are written to a temporary file and renamed, so a frame that
exists is always complete, even while the solver is still running.

| Item | Type | Meaning |
|---|---|---|
| magic | 8 bytes | `EU2DSNAP` |
| version | int32 | format version (1) |
| step | int32 | solver step of this frame |
| time, dx, dy | 3 × float64 | simulation time and cell size |
| nfields | int32 | number of stored fields |
| field ids | nfields × int32 | 0 = ρ, 1 = u, 2 = v, 3 = p |
| nboxes | int32 | number of boxes |
| per box | 4 × int32 | nx, ny, origin_i, origin_j |
| data | float32 | for each box, for each field: nx × ny values, i fastest, bottom row first |

The authoritative description is the comment in `src/io/snapshot.h`; change both together.

## The renderer

`render` turns a directory of frames into a video by piping raw images into ffmpeg (H.264, yuv420p,
`-movflags +faststart` so it plays in browsers). For each video it:

1. **Lists** `frame_*.snap` in the directory, sorted by number.
2. **Pass 1:** reads every frame to find the global minimum and maximum of the field (or of the
   density gradient for schlieren), so the colours mean the same thing in every frame. Skipped
   when `--range` is given.
3. **Composes** each frame into one image: every box is placed at its origin in a global cell grid;
   cells without a box are drawn dark gray. Image row 0 is the top of the domain. The image has the
   physical aspect ratio (W·dx : H·dy) unless `--cells` is given; pixels sample the nearest cell.
4. **Pass 2:** maps values to colours (viridis or gray) and writes the image to ffmpeg.

All frames must have the same box layout; a frame from a different run or resolution stops the
renderer with "Inconsistent frame".

Schlieren images show \(\exp(-k\,|\nabla\rho| / \max|\nabla\rho|)\). Gradients are central
differences on the global grid, so they are continuous across box connections, and one-sided next
to walls and solid regions.

Videos are small because the fields are smooth: a few hundred frames at 432 × 576 take 0.1–0.5 MB
at the default quality. See [Options](options.md#renderer-render) for all flags.

## Reading frames in Python

`tools/snapshot.py` reads a frame into numpy arrays (needs numpy; plotting needs matplotlib):

```python
import sys; sys.path.insert(0, "tools")
import snapshot

s = snapshot.read("results/jet/frame_00100.snap")
print(s["time"], s["step"], s["fields"])    # simulation time, step, stored fields
rho = snapshot.compose(s, "rho")            # whole domain, NaN outside boxes, row 0 = bottom
box = s["boxes"][3]                         # one box: nx, ny, origin_i, origin_j, and its fields
top_rho = box["rho"]                        # shape (ny, nx)
```

From the shell, to plot one frame:

```bash
python tools/snapshot.py results/jet/frame_00100.snap --field rho -o rho.png
```
