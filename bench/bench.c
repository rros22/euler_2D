// Component benchmark: Riemann solver throughput and per-kernel cost in ns per cell.
#define _POSIX_C_SOURCE 199309L
#include "mesh/box.h"
#include "numerics/sweep.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define GAMMA 1.4
#define MIN_TIME 0.5

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

static void bench_riemann(void)
{
    const int n = 1 << 14;
    double *wl = malloc(4 * n * sizeof(double)), *wr = malloc(4 * n * sizeof(double));
    srand(1);
    for (int k = 0; k < n; ++k)
    {
        for (int c = 0; c < 4; ++c)
        {
            double r1 = rand() / (double)RAND_MAX, r2 = rand() / (double)RAND_MAX;
            wl[4 * k + c] = (c == 0 || c == 3) ? 0.5 + r1 : r1 - 0.5;
            wr[4 * k + c] = (c == 0 || c == 3) ? 0.5 + r2 : r2 - 0.5;
        }
    }
    double flux[4], sink = 0.0;
    long calls = 0;
    double start = now();
    while (now() - start < MIN_TIME)
    {
        for (int k = 0; k < n; ++k)
        {
            calculate_flux_1D(flux, &wl[4 * k], &wr[4 * k], GAMMA, 'x');
            sink += flux[0];
        }
        calls += n;
    }
    double el = now() - start;
    printf("riemann  calculate_flux_1D: %8.1f ns/call   (checksum %g)\n", 1e9 * el / calls, sink);
    free(wl);
    free(wr);
}

#define TIME_KERNEL(result, stmt)            \
    do                                       \
    {                                        \
        int reps_ = 0;                       \
        double t0_ = now();                  \
        while (reps_ < 1 || now() - t0_ < MIN_TIME) \
        {                                    \
            stmt;                            \
            ++reps_;                         \
        }                                    \
        result = (now() - t0_) / reps_;      \
    } while (0)

static void reset(struct Box *b)
{
    init_smooth(b);
    box_prim_to_cons(b, GAMMA);
    box_apply_boundary_x_edges(b);
    box_apply_boundary_y_edges(b);
}

static void bench_size(int nx, int ny)
{
    struct Box *b = box_create(nx, ny);
    for (int e = 0; e < 4; ++e)
        box_set_boundary_transmissive(b, e);
    reset(b);

    const struct SchemeConfig first = {RECON_FIRST_ORDER, LIMITER_VAN_LEER};
    const struct SchemeConfig muscl = {RECON_MUSCL_HANCOCK, LIMITER_VAN_LEER};
    double h = 1.0 / nx;
    double dt = 0.4 * h / box_max_wave_speed(b, GAMMA);
    double cells = (double)nx * ny;
    double t_x, t_y, t_mx, t_my, t_c2p, t_bc, t_ws;

    // In-place sweeps with fused cons->prim, as used by solver_step
    TIME_KERNEL(t_x, godunov_x_sweep_cfg(b, b, dt, h, GAMMA, &first, NULL));
    reset(b);
    TIME_KERNEL(t_y, godunov_y_sweep_cfg(b, b, dt, h, GAMMA, &first, NULL));
    reset(b);
    TIME_KERNEL(t_mx, godunov_x_sweep_cfg(b, b, dt, h, GAMMA, &muscl, NULL));
    reset(b);
    TIME_KERNEL(t_my, godunov_y_sweep_cfg(b, b, dt, h, GAMMA, &muscl, NULL));
    reset(b);
    TIME_KERNEL(t_c2p, box_cons_to_prim(b, GAMMA));
    TIME_KERNEL(t_bc, (box_apply_boundary_x_edges(b), box_apply_boundary_y_edges(b)));
    TIME_KERNEL(t_ws, box_max_wave_speed(b, GAMMA));

    printf("%5d x %-5d  first-order x %6.1f y %6.1f   MUSCL-Hancock x %6.1f y %6.1f   "
           "c2p %4.2f  bc %4.2f  wave_speed %4.2f  ns/cell\n",
           nx, ny, 1e9 * t_x / cells, 1e9 * t_y / cells, 1e9 * t_mx / cells, 1e9 * t_my / cells,
           1e9 * t_c2p / cells, 1e9 * t_bc / cells, 1e9 * t_ws / cells);

    box_destroy(b);
}

int main(void)
{
    int sizes[][2] = {{64, 64}, {256, 256}, {1024, 1024}, {4096, 256}};
    bench_riemann();
    for (unsigned s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
        bench_size(sizes[s][0], sizes[s][1]);
    return 0;
}
