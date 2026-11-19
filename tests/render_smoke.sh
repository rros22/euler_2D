#!/bin/sh
# End-to-end smoke test: small jet run -> snapshots -> PPM still (and an mp4 if ffmpeg is available).
# Usage: tests/render_smoke.sh BIN_DIR SCRATCH_DIR
set -e
bin=$1
out=$2/smoke
fail() { echo "FAIL render: $1"; exit 1; }

rm -rf "$out"
"$bin/euler_2D" --case=jet --base=8 --t-end=0.1 --out-dt=0.01 --out-dir="$out" --fields=rho,p > /dev/null
n=$(ls "$out"/frame_*.snap | wc -l | tr -d ' ')
[ "$n" -eq 11 ] || fail "expected 11 frames, got $n"

# base 8: 40 x 64 cells with dx = 1/40, dy = 1/64, so a square 64 x 64 image
"$bin/render" "$out" --still=10 -o "$out/still.ppm" 2> /dev/null || fail "still"
size=$(wc -c < "$out/still.ppm" | tr -d ' ')
[ "$size" -eq $((13 + 64 * 64 * 3)) ] || fail "still.ppm has $size bytes"
# top-left pixel lies outside the boxes: background colour
set -- $(od -An -tu1 -j13 -N3 "$out/still.ppm")
[ "$1 $2 $3" = "40 40 40" ] || fail "background pixel is '$*'"

if command -v ffmpeg > /dev/null 2>&1; then
    "$bin/render" "$out" --schlieren -o "$out/schlieren.mp4" 2> /dev/null || fail "video"
    [ -s "$out/schlieren.mp4" ] || fail "empty video"
    echo "PASS render: frames, still and video"
else
    echo "PASS render: frames and still (SKIP video: ffmpeg not found)"
fi
rm -rf "$out"
