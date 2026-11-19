// Memory-layout microbenchmark: identical first-order kernels (flux once per face, in-place update,
// fused cons->prim) on SoA with separate mallocs, SoA in one padded block, and AoS.
// Run with the exact Riemann flux (what the solver uses) and a cheap Rusanov flux (bandwidth-bound case).
#define _POSIX_C_SOURCE 200112L
#include "mesh/box.h"
#include "numerics/sweep.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define GAMMA 1.4
#define MIN_TIME 0.4

typedef void (*flux_fn)(double *flux, double *l, double *r, double gamma, char dir);

struct Grid
{
    int nx, ny, stride;
    double *w[4], *u[4]; // SoA
    double *aw, *au;     // AoS: 4 doubles per cell
    void *mem[9];
    int nmem;
};

#define GIDX(g, i, j) ((i) + (g)->stride * (j))

#define KNAME(x) soa_##x
#define W(g, c, idx) ((g)->w[c][idx])
#define U(g, c, idx) ((g)->u[c][idx])
#include "bench_layout_kernel.h"
#undef KNAME
#undef W
#undef U

#define KNAME(x) aos_##x
#define W(g, c, idx) ((g)->aw[4 * (idx) + (c)])
#define U(g, c, idx) ((g)->au[4 * (idx) + (c)])
#include "bench_layout_kernel.h"
#undef KNAME
#undef W
#undef U

enum Layout
{
    SOA_SEPARATE,
    SOA_BLOCK,
    AOS
};
static const char *layout_names[] = {"SoA separate", "SoA block", "AoS"};

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void *xalloc(struct Grid *g, size_t bytes)
{
    void *p = NULL;
    if (posix_memalign(&p, 64, bytes) != 0)
    {
        fprintf(stderr, "allocation failed\n");
        exit(1);
    }
    g->mem[g->nmem++] = p;
    return p;
}

static void grid_create(struct Grid *g, int nx, int ny, enum Layout layout)
{
    g->nx = nx;
    g->ny = ny;
    g->stride = nx + 2;
    g->nmem = 0;
    size_t size = (size_t)(nx + 2) * (ny + 2);
    if (layout == SOA_SEPARATE)
    {
        for (int c = 0; c < 4; ++c)
        {
            g->w[c] = xalloc(g, size * sizeof(double));
            g->u[c] = xalloc(g, size * sizeof(double));
        }
    }
    else if (layout == SOA_BLOCK)
    {
        size_t bytes = (size * sizeof(double) + BOX_PAGE - 1) / BOX_PAGE * BOX_PAGE;
        size_t stride = bytes + BOX_ARRAY_SKEW;
        char *base = xalloc(g, 8 * stride);
        for (int c = 0; c < 4; ++c)
        {
            g->w[c] = (double *)(base + c * stride);
            g->u[c] = (double *)(base + (4 + c) * stride);
        }
    }
    else
    {
        g->aw = xalloc(g, 4 * size * sizeof(double));
        g->au = xalloc(g, 4 * size * sizeof(double));
    }

    const double pi = 3.14159265358979323846;
    for (int j = 0; j <= ny + 1; ++j)
        for (int i = 0; i <= nx + 1; ++i)
        {
            int ic = i < 1 ? 1 : (i > nx ? nx : i), jc = j < 1 ? 1 : (j > ny ? ny : j);
            double x = (ic - 0.5) / nx, y = (jc - 0.5) / ny;
            double w[4] = {1.0 + 0.3 * sin(2 * pi * x) * cos(2 * pi * y), 0.5 * sin(4 * pi * y),
                           0.5 * cos(4 * pi * x), 1.0 + 0.3 * cos(2 * pi * (x + y))};
            double u[4] = {w[0], w[0] * w[1], w[0] * w[2],
                           0.5 * w[0] * (w[1] * w[1] + w[2] * w[2]) + w[3] / (GAMMA - 1.0)};
            int idx = GIDX(g, i, j);
            for (int c = 0; c < 4; ++c)
            {
                if (layout == AOS)
                {
                    g->aw[4 * idx + c] = w[c];
                    g->au[4 * idx + c] = u[c];
                }
                else
                {
                    g->w[c][idx] = w[c];
                    g->u[c][idx] = u[c];
                }
            }
        }
}

static void grid_destroy(struct Grid *g)
{
    for (int k = 0; k < g->nmem; ++k)
        free(g->mem[k]);
}

static void rusanov(double *flux, double *l, double *r, double gamma, char dir)
{
    int n = (dir == 'x') ? 1 : 2;
    double fl[4], fr[4], ul[4], ur[4], smax = 0.0;
    double *w[2] = {l, r}, *f[2] = {fl, fr}, *u[2] = {ul, ur};
    for (int s = 0; s < 2; ++s)
    {
        double rho = w[s][0], vn = w[s][n], p = w[s][3];
        double E = 0.5 * rho * (w[s][1] * w[s][1] + w[s][2] * w[s][2]) + p / (gamma - 1.0);
        u[s][0] = rho;
        u[s][1] = rho * w[s][1];
        u[s][2] = rho * w[s][2];
        u[s][3] = E;
        f[s][0] = rho * vn;
        f[s][1] = rho * vn * w[s][1];
        f[s][2] = rho * vn * w[s][2];
        f[s][n] += p;
        f[s][3] = vn * (E + p);
        double sp = fabs(vn) + sqrt(gamma * p / rho);
        if (sp > smax)
            smax = sp;
    }
    for (int c = 0; c < 4; ++c)
        flux[c] = 0.5 * (fl[c] + fr[c]) - 0.5 * smax * (ur[c] - ul[c]);
}

static void exact(double *flux, double *l, double *r, double gamma, char dir)
{
    calculate_flux_1D(flux, l, r, gamma, dir);
}

// kernel: 0 = x-sweep, 1 = y-sweep, 2 = cons->prim only
static double time_kernel(int nx, int ny, enum Layout layout, int kernel, flux_fn flux)
{
    struct Grid g;
    grid_create(&g, nx, ny, layout);
    double *F = malloc(8 * (size_t)(nx + 1) * sizeof(double));
    double dt_h = 0.05;
    int reps = 0;
    double t0 = now();
    while (reps < 2 || now() - t0 < MIN_TIME)
    {
        if (layout == AOS)
        {
            if (kernel == 0)
                aos_x_sweep(&g, dt_h, GAMMA, flux, F);
            else if (kernel == 1)
                aos_y_sweep(&g, dt_h, GAMMA, flux, F);
            else
                aos_convert_all(&g, GAMMA);
        }
        else
        {
            if (kernel == 0)
                soa_x_sweep(&g, dt_h, GAMMA, flux, F);
            else if (kernel == 1)
                soa_y_sweep(&g, dt_h, GAMMA, flux, F);
            else
                soa_convert_all(&g, GAMMA);
        }
        ++reps;
    }
    double t = (now() - t0) / reps;
    free(F);
    grid_destroy(&g);
    return 1e9 * t / ((double)nx * ny);
}

int main(void)
{
    int sizes[][2] = {{256, 256}, {1024, 1024}, {2048, 2048}, {4096, 256}};
    struct
    {
        const char *name;
        flux_fn fn;
    } fluxes[] = {{"exact", exact}, {"rusanov", rusanov}};

    printf("ns per cell; kernels: x/y = sweep with flux once per face + update + cons->prim, c2p = cons->prim only\n");
    printf("%-11s %-8s %-13s %9s %9s %9s\n", "size", "flux", "layout", "x-sweep", "y-sweep", "c2p");
    for (unsigned s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
        for (int f = 0; f < 2; ++f)
            for (int l = 0; l < 3; ++l)
            {
                int nx = sizes[s][0], ny = sizes[s][1];
                double tx = time_kernel(nx, ny, l, 0, fluxes[f].fn);
                double ty = time_kernel(nx, ny, l, 1, fluxes[f].fn);
                double tc = f == 0 ? time_kernel(nx, ny, l, 2, NULL) : NAN;
                char size[32];
                snprintf(size, sizeof(size), "%dx%d", nx, ny);
                printf("%-11s %-8s %-13s %9.2f %9.2f %9.2f\n", size, fluxes[f].name, layout_names[l], tx, ty, tc);
            }
    return 0;
}
