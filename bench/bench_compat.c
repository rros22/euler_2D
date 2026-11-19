// Sweep benchmark that uses only the original public API (box_*, godunov_x/y_sweep, box_swap):
// the old-style loop sweeps into a temp box, converts, and swaps. Against the original sources
// (git history) it compiles with the include paths box/box.h and godunov/godunov.h.
#define _POSIX_C_SOURCE 199309L
#include "mesh/box.h"
#include "numerics/sweep.h"
#include <math.h>
#include <stdio.h>
#include <time.h>

#define GAMMA 1.4

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void init_smooth(struct Box *b)
{
    const double pi = 3.14159265358979323846;
    for (int j = 1; j <= b->ny; ++j)
        for (int i = 1; i <= b->nx; ++i)
        {
            double x = (i - 0.5) / b->nx, y = (j - 0.5) / b->ny;
            int idx = INDEX(i, j, b->nx);
            b->w_1[idx] = 1.0 + 0.3 * sin(2 * pi * x) * cos(2 * pi * y);
            b->w_2[idx] = 0.5 * sin(4 * pi * y);
            b->w_3[idx] = 0.5 * cos(4 * pi * x);
            b->w_4[idx] = 1.0 + 0.3 * cos(2 * pi * (x + y));
        }
}

static void pass(struct Box **curr, struct Box **temp, int dir, double dt, double h)
{
    if (dir == 0)
    {
        box_apply_boundary_x_edges(*curr);
        godunov_x_sweep(*curr, *temp, dt, h, GAMMA);
    }
    else
    {
        box_apply_boundary_y_edges(*curr);
        godunov_y_sweep(*curr, *temp, dt, h, GAMMA);
    }
    box_cons_to_prim(*temp, GAMMA);
    box_swap(curr, temp);
}

static void run(int nx, int ny)
{
    struct Box *b = box_create(nx, ny), *t = box_create(nx, ny);
    for (int e = 0; e < 4; ++e)
        box_set_boundary_transmissive(b, e);
    box_copy_boundaries(t, b);
    init_smooth(b);
    box_prim_to_cons(b, GAMMA);

    double h = 1.0 / nx;
    double dt = 0.4 * h / box_max_wave_speed(b, GAMMA);
    double cells = (double)nx * ny;

    double tx = 0.0, ty = 0.0;
    int reps = 0;
    double start = now();
    while (reps < 1 || now() - start < 1.0)
    {
        double t0 = now();
        pass(&b, &t, 0, dt, h);
        double t1 = now();
        pass(&b, &t, 1, dt, h);
        double t2 = now();
        tx += t1 - t0;
        ty += t2 - t1;
        ++reps;
    }
    printf("%5d x %-5d  reps %4d   x-pass %8.1f ns/cell   y-pass %8.1f ns/cell\n",
           nx, ny, reps, 1e9 * tx / (reps * cells), 1e9 * ty / (reps * cells));

    box_destroy(b);
    box_destroy(t);
}

int main(void)
{
    int sizes[][2] = {{64, 64}, {256, 256}, {1024, 1024}, {4096, 256}};
    printf("bench_compat: old-API pass = boundaries + sweep + cons_to_prim + swap\n");
    for (unsigned s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
        run(sizes[s][0], sizes[s][1]);
    return 0;
}
