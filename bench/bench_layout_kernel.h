// Kernel body instantiated once per layout by bench_layout.c.
// Requires: KNAME(x) name mangling, W(g, c, idx) and U(g, c, idx) lvalue accessors.

static void KNAME(convert_row)(struct Grid *g, int j, double gamma)
{
    for (int i = 1; i <= g->nx; ++i)
    {
        int idx = GIDX(g, i, j);
        double w[4];
        cons_to_prim_cell(U(g, 0, idx), U(g, 1, idx), U(g, 2, idx), U(g, 3, idx), gamma, w);
        W(g, 0, idx) = w[0];
        W(g, 1, idx) = w[1];
        W(g, 2, idx) = w[2];
        W(g, 3, idx) = w[3];
    }
}

static void KNAME(x_sweep)(struct Grid *g, double dt_h, double gamma, flux_fn flux, double *F)
{
    int nx = g->nx, n = nx + 1;
    for (int j = 1; j <= g->ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            int a = GIDX(g, i, j), b = GIDX(g, i + 1, j);
            double l[4] = {W(g, 0, a), W(g, 1, a), W(g, 2, a), W(g, 3, a)};
            double r[4] = {W(g, 0, b), W(g, 1, b), W(g, 2, b), W(g, 3, b)};
            double f[4];
            flux(f, l, r, gamma, 'x');
            for (int c = 0; c < 4; ++c)
                F[c * n + i] = f[c];
        }
        for (int i = 1; i <= nx; ++i)
        {
            int idx = GIDX(g, i, j);
            for (int c = 0; c < 4; ++c)
                U(g, c, idx) = U(g, c, idx) + dt_h * (F[c * n + i - 1] - F[c * n + i]);
        }
        KNAME(convert_row)(g, j, gamma);
    }
}

static void KNAME(y_sweep)(struct Grid *g, double dt_h, double gamma, flux_fn flux, double *F)
{
    int nx = g->nx;
    double *Fm = F, *Fp = F + 4 * (nx + 1);
    for (int j = 0; j <= g->ny; ++j)
    {
        for (int i = 1; i <= nx; ++i)
        {
            int a = GIDX(g, i, j), b = GIDX(g, i, j + 1);
            double l[4] = {W(g, 0, a), W(g, 1, a), W(g, 2, a), W(g, 3, a)};
            double r[4] = {W(g, 0, b), W(g, 1, b), W(g, 2, b), W(g, 3, b)};
            double f[4];
            flux(f, l, r, gamma, 'y');
            for (int c = 0; c < 4; ++c)
                Fp[c * nx + i - 1] = f[c];
        }
        if (j >= 1)
        {
            for (int i = 1; i <= nx; ++i)
            {
                int idx = GIDX(g, i, j);
                for (int c = 0; c < 4; ++c)
                    U(g, c, idx) = U(g, c, idx) + dt_h * (Fm[c * nx + i - 1] - Fp[c * nx + i - 1]);
            }
            KNAME(convert_row)(g, j, gamma);
        }
        double *t = Fm;
        Fm = Fp;
        Fp = t;
    }
}

static void KNAME(convert_all)(struct Grid *g, double gamma)
{
    for (int j = 1; j <= g->ny; ++j)
        KNAME(convert_row)(g, j, gamma);
}
