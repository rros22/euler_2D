// Physics tests: ghost layers, conservation, convergence order, Sod/Toro shock tubes, x/y symmetry.
#include "test_util.h"
#include "numerics/riemann_exact.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GAMMA 1.4
#define PI 3.14159265358979323846

static int failures = 0;

static void check(int ok, const char *name, const char *fmt, double value)
{
    printf("%s %-52s ", ok ? "PASS" : "FAIL", name);
    printf(fmt, value);
    printf("\n");
    if (!ok)
        ++failures;
}

struct Scheme
{
    const char *name;
    struct SchemeConfig cfg;
};

static const struct Scheme schemes[] = {
    {"first-order", {RECON_FIRST_ORDER, LIMITER_MINMOD}},
    {"MUSCL-Hancock minmod", {RECON_MUSCL_HANCOCK, LIMITER_MINMOD}},
    {"MUSCL-Hancock van Leer", {RECON_MUSCL_HANCOCK, LIMITER_VAN_LEER}},
    {"MUSCL-Hancock MC", {RECON_MUSCL_HANCOCK, LIMITER_MC}},
};
#define NUM_SCHEMES ((int)(sizeof(schemes) / sizeof(schemes[0])))

// ---------------------------------------------------------------------------------------------
// Every ghost layer of every boundary type holds the expected source value

static void test_ghost_layers(void)
{
    const int nx = 5, ny = 4;
    struct Box *b = box_create(nx, ny), *nb = box_create(nx, ny);
    for (int j = 1; j <= ny; ++j)
        for (int i = 1; i <= nx; ++i)
        {
            int idx = INDEX(i, j, nx);
            b->w_1[idx] = 1.0 + i + 10.0 * j;
            b->w_2[idx] = 0.1 * i;
            b->w_3[idx] = 0.2 * j;
            b->w_4[idx] = 2.0 + i * j;
            nb->w_1[idx] = 100.0 + i + 10.0 * j;
            nb->w_2[idx] = nb->w_3[idx] = nb->w_4[idx] = 1.0;
        }
    double inflow[4] = {7.0, 8.0, 9.0, 10.0};

    box_set_boundary_reflective(b, EDGE_LEFT);
    box_set_boundary_connected(b, EDGE_RIGHT, nb, EDGE_LEFT);
    box_set_boundary_transmissive(b, EDGE_BOTTOM);
    box_set_boundary_inflow(b, inflow, EDGE_TOP);
    box_apply_boundary_x_edges(b);
    box_apply_boundary_y_edges(b);

    int bad = 0;
    for (int g = 1; g <= BOX_NG; ++g)
    {
        for (int j = 1; j <= ny; ++j)
        {
            int gl = INDEX(1 - g, j, nx), ml = INDEX(g, j, nx);
            bad += b->w_1[gl] != b->w_1[ml] || b->w_2[gl] != -b->w_2[ml] || b->w_3[gl] != b->w_3[ml];
            int gr = INDEX(nx + g, j, nx), sr = INDEX(g, j, nx);
            bad += b->w_1[gr] != nb->w_1[sr];
        }
        for (int i = 1; i <= nx; ++i)
        {
            int gb = INDEX(i, 1 - g, nx), eb = INDEX(i, 1, nx);
            bad += b->w_1[gb] != b->w_1[eb] || b->w_3[gb] != b->w_3[eb];
            int gt = INDEX(i, ny + g, nx);
            bad += b->w_1[gt] != 7.0 || b->w_4[gt] != 10.0;
        }
    }
    check(bad == 0, "ghost layers: reflective/connected/transmissive/inflow", "mismatches %.0f", (double)bad);

    box_set_boundary_reflective(b, EDGE_TOP);
    box_apply_boundary_y_edges(b);
    bad = 0;
    for (int g = 1; g <= BOX_NG; ++g)
        for (int i = 1; i <= nx; ++i)
        {
            int gt = INDEX(i, ny + g, nx), mt = INDEX(i, ny + 1 - g, nx);
            bad += b->w_1[gt] != b->w_1[mt] || b->w_3[gt] != -b->w_3[mt] || b->w_2[gt] != b->w_2[mt];
        }
    check(bad == 0, "ghost layers: reflective top flips v", "mismatches %.0f", (double)bad);

    box_destroy(b);
    box_destroy(nb);
}

// ---------------------------------------------------------------------------------------------
// Closed 4-box domain (connected interior edges, reflective outer edges): mass and energy conserved

static void test_conservation(void)
{
    for (int s = 0; s < NUM_SCHEMES; ++s)
    {
        struct Setup st;
        setup_four_box(&st, 10, 1, GAMMA);
        double m0 = setup_total(&st, 1), e0 = setup_total(&st, 4);
        setup_run(&st, GAMMA, 0.4, 1e30, 100, &schemes[s].cfg);
        double dm = fabs(setup_total(&st, 1) - m0) / m0;
        double de = fabs(setup_total(&st, 4) - e0) / e0;
        char name[96];
        snprintf(name, sizeof(name), "conservation closed 4-box, %s: mass", schemes[s].name);
        check(dm < 1e-13, name, "rel. error %.3e", dm);
        snprintf(name, sizeof(name), "conservation closed 4-box, %s: energy", schemes[s].name);
        check(de < 1e-13, name, "rel. error %.3e", de);
        check(setup_is_physical(&st), "  positivity", "%.0f", 0.0);
        setup_destroy(&st);
    }
}

// ---------------------------------------------------------------------------------------------
// Quasi-1D problems on a single box: N cells along dir (0 = x, 1 = y), 4 cells across

typedef void (*init_fn)(double a, double b, const void *ctx, double *w); // cell [a, b] -> rho, u_n, u_t, p

static void setup_line(struct Setup *s, int N, int dir, int periodic, init_fn f, const void *ctx)
{
    int nx = dir == 0 ? N : 4, ny = dir == 0 ? 4 : N;
    double h = 1.0 / N;
    setup_single_box(s, nx, ny, h, h);
    struct Box *b = s->boxes[0];
    if (periodic)
    {
        box_set_boundary_connected(b, EDGE_LEFT, b, EDGE_RIGHT);
        box_set_boundary_connected(b, EDGE_RIGHT, b, EDGE_LEFT);
        box_set_boundary_connected(b, EDGE_BOTTOM, b, EDGE_TOP);
        box_set_boundary_connected(b, EDGE_TOP, b, EDGE_BOTTOM);
    }
    else
    {
        for (int e = 0; e < 4; ++e)
            box_set_boundary_transmissive(b, e);
    }
    for (int j = 1; j <= ny; ++j)
        for (int i = 1; i <= nx; ++i)
        {
            int k = dir == 0 ? i : j;
            double w[4];
            f((k - 1) * h, k * h, ctx, w);
            int idx = INDEX(i, j, nx);
            b->w_1[idx] = w[0];
            b->w_2[idx] = dir == 0 ? w[1] : w[2];
            b->w_3[idx] = dir == 0 ? w[2] : w[1];
            b->w_4[idx] = w[3];
        }
    setup_finish(s, GAMMA);
    // Make w exactly what every sweep produces from u, so x and y runs start from identical states
    // regardless of which direction's (trivial) sweep happens to run first
    box_cons_to_prim(b, GAMMA);
}

// Density along the line (first transverse row/column)
static void line_density(const struct Setup *s, int N, int dir, double *rho)
{
    const struct Box *b = s->boxes[0];
    for (int k = 1; k <= N; ++k)
        rho[k - 1] = b->w_1[dir == 0 ? INDEX(k, 1, b->nx) : INDEX(1, k, b->nx)];
}

// --- smooth density wave advected with u = 1, p = 1; exact cell averages

static void init_wave(double a, double b, const void *ctx, double *w)
{
    (void)ctx;
    w[0] = 1.0 + 0.2 * (cos(2 * PI * a) - cos(2 * PI * b)) / (2 * PI * (b - a));
    w[1] = 1.0;
    w[2] = 0.0;
    w[3] = 1.0;
}

static double advection_error(int N, int dir, const struct SchemeConfig *cfg, double cfl)
{
    struct Setup s;
    setup_line(&s, N, dir, 1, init_wave, NULL);
    setup_run(&s, GAMMA, cfl, 1.0, 1000000, cfg); // one period: exact solution = initial data
    double *rho = malloc(N * sizeof(double)), err = 0.0;
    line_density(&s, N, dir, rho);
    for (int k = 1; k <= N; ++k)
    {
        double w[4];
        init_wave((k - 1.0) / N, (double)k / N, NULL, w);
        err += fabs(rho[k - 1] - w[0]) / N;
    }
    free(rho);
    setup_destroy(&s);
    return err;
}

static void test_convergence(void)
{
    const int Ns[] = {32, 64, 128, 256};
    const int nN = 4;
    // Minimum observed L1 order on the finest pair
    const double min_order[NUM_SCHEMES] = {0.85, 1.6, 1.8, 1.8};

    printf("\nSmooth density wave, L1 error after one period (CFL 0.8)\n");
    printf("%-24s %4s %12s %12s %7s\n", "scheme", "N", "L1 (x)", "L1 (y)", "order");
    for (int s = 0; s < NUM_SCHEMES; ++s)
    {
        double ex[4], ey[4], order = 0.0;
        int symmetric = 1;
        for (int n = 0; n < nN; ++n)
        {
            ex[n] = advection_error(Ns[n], 0, &schemes[s].cfg, 0.8);
            ey[n] = advection_error(Ns[n], 1, &schemes[s].cfg, 0.8);
            symmetric &= (ex[n] == ey[n]);
            if (n > 0)
                order = log(ex[n - 1] / ex[n]) / log(2.0);
            if (n == 0)
                printf("%-24s %4d %12.4e %12.4e %7s\n", schemes[s].name, Ns[n], ex[n], ey[n], "-");
            else
                printf("%-24s %4d %12.4e %12.4e %7.2f\n", "", Ns[n], ex[n], ey[n], order);
        }
        char name[96];
        snprintf(name, sizeof(name), "advection order, %s", schemes[s].name);
        check(order > min_order[s], name, "order %.2f", order);
        snprintf(name, sizeof(name), "advection x/y errors identical, %s", schemes[s].name);
        check(symmetric, name, "%.0f", 0.0);
    }
    printf("\n");
}

// --- Riemann problems (Toro, Riemann Solvers and Numerical Methods, Table 4.1)

struct RiemannTest
{
    const char *name;
    double wl[3], wr[3]; // rho, u, p
    double x0, t_end;
};

static const struct RiemannTest riemann_tests[] = {
    {"Sod", {1.0, 0.0, 1.0}, {0.125, 0.0, 0.1}, 0.5, 0.2},
    {"Toro 1", {1.0, 0.75, 1.0}, {0.125, 0.0, 0.1}, 0.3, 0.2},
    {"Toro 2 (123)", {1.0, -2.0, 0.4}, {1.0, 2.0, 0.4}, 0.5, 0.15},
    {"Toro 3", {1.0, 0.0, 1000.0}, {1.0, 0.0, 0.01}, 0.5, 0.012},
    {"Toro 4", {5.99924, 19.5975, 460.894}, {5.99242, -6.19633, 46.0950}, 0.4, 0.035},
    {"Toro 5", {1.0, -19.59745, 1000.0}, {1.0, -19.59745, 0.01}, 0.8, 0.012},
};
#define NUM_RIEMANN ((int)(sizeof(riemann_tests) / sizeof(riemann_tests[0])))

static void init_riemann(double a, double b, const void *ctx, double *w)
{
    const struct RiemannTest *rt = ctx;
    const double *s = (0.5 * (a + b) < rt->x0) ? rt->wl : rt->wr;
    w[0] = s[0];
    w[1] = s[1];
    w[2] = 0.0;
    w[3] = s[2];
}

static double exact_density(const struct RiemannTest *rt, double x)
{
    double wl[3] = {rt->wl[0], rt->wl[1], rt->wl[2]}, wr[3] = {rt->wr[0], rt->wr[1], rt->wr[2]};
    double sol[3];
    riemann_solver(sol, wl, wr, (x - rt->x0) / rt->t_end, GAMMA, 1e-12);
    return sol[0];
}

// Returns L1 density error; physical = 0 if any cell has non-finite or non-positive rho or p
static double riemann_run(const struct RiemannTest *rt, int N, int dir, const struct SchemeConfig *cfg,
                          double cfl, int *physical, double *rho_out)
{
    struct Setup s;
    setup_line(&s, N, dir, 0, init_riemann, rt);
    setup_run(&s, GAMMA, cfl, rt->t_end, 1000000, cfg);
    *physical = setup_is_physical(&s);
    double *rho = rho_out ? rho_out : malloc(N * sizeof(double)), err = 0.0;
    line_density(&s, N, dir, rho);
    for (int k = 1; k <= N; ++k)
        err += fabs(rho[k - 1] - exact_density(rt, (k - 0.5) / N)) / N;
    if (!rho_out)
        free(rho);
    setup_destroy(&s);
    return err;
}

static void test_riemann_problems(void)
{
    const int N = 200;
    const double cfls[] = {0.4, 0.9};

    printf("Riemann problems, N = %d, L1 density error\n", N);
    printf("%-14s %-24s %10s %10s\n", "test", "scheme", "CFL 0.4", "CFL 0.9");
    for (int r = 0; r < NUM_RIEMANN; ++r)
    {
        double err_fo = 0.0;
        for (int s = 0; s < NUM_SCHEMES; ++s)
        {
            double err[2];
            int all_physical = 1;
            for (int c = 0; c < 2; ++c)
            {
                int physical;
                err[c] = riemann_run(&riemann_tests[r], N, 0, &schemes[s].cfg, cfls[c], &physical, NULL);
                all_physical &= physical && isfinite(err[c]);
            }
            printf("%-14s %-24s %10.4e %10.4e\n", s == 0 ? riemann_tests[r].name : "", schemes[s].name, err[0], err[1]);
            char name[96];
            snprintf(name, sizeof(name), "  %s, %s: finite, rho > 0, p > 0", riemann_tests[r].name, schemes[s].name);
            if (!all_physical)
                check(0, name, "%.0f", 0.0);
            if (s == 0)
                err_fo = err[0];
            else if (r == 0)
            {
                snprintf(name, sizeof(name), "  Sod: %s error below first order", schemes[s].name);
                check(err[0] < err_fo, name, "ratio %.3f", err[0] / err_fo);
            }
        }
    }
    printf("\n");
}

// Sod along x and along y give bitwise identical density profiles
static void test_xy_symmetry(void)
{
    const int N = 100;
    for (int s = 0; s < NUM_SCHEMES; ++s)
    {
        double rx[N], ry[N];
        int px, py;
        riemann_run(&riemann_tests[0], N, 0, &schemes[s].cfg, 0.8, &px, rx);
        riemann_run(&riemann_tests[0], N, 1, &schemes[s].cfg, 0.8, &py, ry);
        double max_diff = 0.0;
        for (int k = 0; k < N; ++k)
            max_diff = fmax(max_diff, fabs(rx[k] - ry[k]));
        char name[96];
        snprintf(name, sizeof(name), "Sod x vs y transposed identical, %s", schemes[s].name);
        check(max_diff == 0.0 && px && py, name, "max |diff| %.3e", max_diff);
    }
}

int main(void)
{
    test_ghost_layers();
    test_conservation();
    test_xy_symmetry();
    test_convergence();
    test_riemann_problems();
    printf("%s: %d failure(s)\n", failures ? "FAILED" : "OK", failures);
    return failures ? 1 : 0;
}
