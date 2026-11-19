"""Reader for euler_2D binary snapshots (frame_*.snap), see src/io/snapshot.h for the format.

    import snapshot
    s = snapshot.read("results/jet/frame_00100.snap")
    rho = snapshot.compose(s, "rho")      # global 2D array, NaN outside the boxes, row 0 = bottom
    s["boxes"][3]["rho"]                  # one box, shape (ny, nx)

Command line: python tools/snapshot.py FRAME [--field rho] [-o out.png]
"""
import sys

import numpy as np

FIELD_NAMES = ["rho", "u", "v", "p"]
MAGIC = b"EU2DSNAP"


def read(path):
    """Return a dict with step, time, dx, dy, fields and boxes (each with nx, ny, origin_i, origin_j
    and one (ny, nx) float32 array per field)."""
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:8] != MAGIC:
        raise ValueError(f"{path}: not a snapshot")
    pos = 8
    version, step = np.frombuffer(raw, "<i4", 2, pos)
    pos += 8
    if version != 1:
        raise ValueError(f"{path}: unsupported version {version}")
    time, dx, dy = np.frombuffer(raw, "<f8", 3, pos)
    pos += 24
    nfields = int(np.frombuffer(raw, "<i4", 1, pos)[0])
    pos += 4
    fields = [FIELD_NAMES[k] for k in np.frombuffer(raw, "<i4", nfields, pos)]
    pos += 4 * nfields
    nboxes = int(np.frombuffer(raw, "<i4", 1, pos)[0])
    pos += 4
    geom = np.frombuffer(raw, "<i4", 4 * nboxes, pos).reshape(nboxes, 4)
    pos += 16 * nboxes

    boxes = []
    for nx, ny, oi, oj in geom.tolist():
        box = {"nx": nx, "ny": ny, "origin_i": oi, "origin_j": oj}
        for name in fields:
            box[name] = np.frombuffer(raw, "<f4", nx * ny, pos).reshape(ny, nx)
            pos += 4 * nx * ny
        boxes.append(box)
    return {"step": int(step), "time": float(time), "dx": float(dx), "dy": float(dy),
            "fields": fields, "boxes": boxes}


def compose(snap, field="rho"):
    """All boxes placed in one global array (ny_total, nx_total); NaN where there is no box."""
    boxes = snap["boxes"]
    W = max(b["origin_i"] + b["nx"] for b in boxes)
    H = max(b["origin_j"] + b["ny"] for b in boxes)
    out = np.full((H, W), np.nan, dtype=np.float32)
    for b in boxes:
        out[b["origin_j"]:b["origin_j"] + b["ny"], b["origin_i"]:b["origin_i"] + b["nx"]] = b[field]
    return out


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("frame")
    ap.add_argument("--field", default="rho")
    ap.add_argument("-o", "--output", help="save the figure instead of showing it")
    args = ap.parse_args()

    s = read(args.frame)
    img = compose(s, args.field)
    print(f"step {s['step']}  t = {s['time']:.5f}  {args.field} in [{np.nanmin(img):.6g}, {np.nanmax(img):.6g}]")

    import matplotlib
    if args.output:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    H, W = img.shape
    plt.imshow(img, origin="lower", extent=(0, W * s["dx"], 0, H * s["dy"]), cmap="viridis")
    plt.colorbar(label=args.field)
    plt.title(f"t = {s['time']:.4f}")
    if args.output:
        plt.savefig(args.output, dpi=150)
    else:
        plt.show()
    sys.exit(0)
