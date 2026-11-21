// Render snapshot frames (frame_*.snap) to an H.264 video by piping raw RGB frames into ffmpeg.
// Pass 1 finds a global colour range over all frames, pass 2 composes the boxes into one image per frame.
// Run "render --help" for the options; a single frame can also be written as a PPM image (--still).
// Documented in docs/options.md and docs/output-and-rendering.md; keep usage() and the docs in sync.
#define _POSIX_C_SOURCE 200809L
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io/snapshot.h"

struct Options
{
    const char *dir; // directory with the frames
    char out[1024];  // output file
    const char *ffmpeg;
    const char *preset; // x264 speed/size trade-off
    int field;          // SNAP_FIELD_* to show
    int schlieren;
    double schlieren_k; // schlieren contrast
    int log_scale;
    int gray; // gray instead of viridis colour map (-1 = automatic)
    int have_range;
    double lo, hi; // fixed colour range if have_range
    double scale;  // output height in pixels per cell
    int physical; // pixel aspect from dx, dy (default) instead of one pixel per cell
    int fps, crf, stride;
    int still;          // >= 0: write this frame as a PPM image instead of a video
    unsigned char bg[3]; // colour outside the boxes
};

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s DIR [options]   (DIR contains frame_*.snap)\n"
            "  --field=rho|u|v|p     field to show (default rho)\n"
            "  --schlieren           numerical schlieren exp(-k |grad rho| / max), gray\n"
            "  --schlieren-k=K       schlieren contrast (default 15)\n"
            "  --log                 log10 of the field\n"
            "  --range=LO,HI         fixed colour range (default: global min/max over all frames)\n"
            "  --cmap=viridis|gray   colour map (default viridis, gray for schlieren)\n"
            "  --scale=S             output height = S * cells in y (default 1)\n"
            "  --cells               one pixel per cell (ignore physical dx/dy aspect)\n"
            "  --stride=K            use every K-th frame (default 1)\n"
            "  --fps=N               frame rate (default 30)\n"
            "  --crf=N               x264 quality, larger = smaller file (default 26)\n"
            "  --preset=NAME         x264 preset (default slow)\n"
            "  --ffmpeg=PATH         ffmpeg executable (default ffmpeg)\n"
            "  --still=N             write frame N as a PPM image to -o instead of a video\n"
            "  -o FILE               output (default DIR/<field>.mp4)\n"
            "  --help                show this message\n",
            prog);
}

static int parse_args(int argc, char **argv, struct Options *opt)
{
    for (int a = 1; a < argc; ++a)
    {
        const char *arg = argv[a];
        const char *val = strchr(arg, '=');
        val = val ? val + 1 : "";
        if (strncmp(arg, "--field=", 8) == 0)
        {
            int f[SNAP_MAX_FIELDS];
            if (snapshot_parse_fields(val, f) != 1)
            {
                fprintf(stderr, "Invalid field '%s'\n", val);
                return -1;
            }
            opt->field = f[0];
        }
        else if (strcmp(arg, "--schlieren") == 0)
            opt->schlieren = 1;
        else if (strncmp(arg, "--schlieren-k=", 14) == 0)
            opt->schlieren_k = atof(val);
        else if (strcmp(arg, "--log") == 0)
            opt->log_scale = 1;
        else if (strncmp(arg, "--range=", 8) == 0)
        {
            if (sscanf(val, "%lf,%lf", &opt->lo, &opt->hi) != 2 || !(opt->hi > opt->lo))
            {
                fprintf(stderr, "Invalid range '%s'\n", val);
                return -1;
            }
            opt->have_range = 1;
        }
        else if (strcmp(arg, "--cmap=viridis") == 0)
            opt->gray = 0;
        else if (strcmp(arg, "--cmap=gray") == 0)
            opt->gray = 1;
        else if (strncmp(arg, "--scale=", 8) == 0)
            opt->scale = atof(val);
        else if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0)
            return -1; // prints the usage
        else if (strcmp(arg, "--cells") == 0)
            opt->physical = 0;
        else if (strncmp(arg, "--stride=", 9) == 0)
            opt->stride = atoi(val);
        else if (strncmp(arg, "--fps=", 6) == 0)
            opt->fps = atoi(val);
        else if (strncmp(arg, "--crf=", 6) == 0)
            opt->crf = atoi(val);
        else if (strncmp(arg, "--preset=", 9) == 0)
            opt->preset = val;
        else if (strncmp(arg, "--ffmpeg=", 9) == 0)
            opt->ffmpeg = val;
        else if (strncmp(arg, "--still=", 8) == 0)
            opt->still = atoi(val);
        else if (strcmp(arg, "-o") == 0 && a + 1 < argc)
            snprintf(opt->out, sizeof(opt->out), "%s", argv[++a]);
        else if (arg[0] != '-' && !opt->dir)
            opt->dir = arg;
        else
        {
            fprintf(stderr, "Unknown option %s\n", arg);
            return -1;
        }
    }
    if (!opt->dir || opt->scale <= 0.0 || opt->stride < 1 || opt->fps < 1)
        return -1;
    return 0;
}

static int compare_names(const void *a, const void *b)
{
    return strcmp(*(char *const *)a, *(char *const *)b);
}

// Sorted list of DIR/frame_*.snap paths
static char **list_frames(const char *dir, int *count)
{
    DIR *d = opendir(dir);
    if (!d)
    {
        perror(dir);
        return NULL;
    }
    char **names = NULL;
    int n = 0, cap = 0;
    struct dirent *e;
    while ((e = readdir(d)))
    {
        size_t len = strlen(e->d_name);
        if (strncmp(e->d_name, "frame_", 6) != 0 || len < 11 || strcmp(e->d_name + len - 5, ".snap") != 0)
            continue;
        if (n == cap)
        {
            cap = cap ? 2 * cap : 256;
            char **tmp = realloc(names, cap * sizeof(*names));
            if (!tmp)
                break;
            names = tmp;
        }
        size_t size = strlen(dir) + len + 2;
        names[n] = malloc(size);
        if (!names[n])
            break;
        snprintf(names[n], size, "%s/%s", dir, e->d_name);
        ++n;
    }
    closedir(d);
    qsort(names, n, sizeof(*names), compare_names);
    *count = n;
    return names;
}

// Viridis colour map from a polynomial fit, tabulated at 256 levels
static void build_lut(unsigned char lut[256][3], int gray)
{
    static const double c[7][3] = {
        {0.2777273272234177, 0.005407344544966578, 0.3340998053353061},
        {0.1050930431085774, 1.404613529898575, 1.384590162594685},
        {-0.3308618287255563, 0.214847559468213, 0.09509516302823659},
        {-4.634230498983486, -5.799100973351585, -19.33244095627987},
        {6.228269936347081, 14.17993336680509, 56.69055260068105},
        {4.776384997670288, -13.74514537774601, -65.35303263337234},
        {-5.435455855934631, 4.645852612178535, 26.3124352495832},
    };
    for (int k = 0; k < 256; ++k)
    {
        double t = k / 255.0;
        for (int ch = 0; ch < 3; ++ch)
        {
            double v = t;
            if (!gray)
            {
                v = c[6][ch];
                for (int p = 5; p >= 0; --p)
                    v = c[p][ch] + t * v;
            }
            v = v < 0.0 ? 0.0 : v > 1.0 ? 1.0 : v;
            lut[k][ch] = (unsigned char)(255.0 * v + 0.5);
        }
    }
}

// Global cell grid holding the displayed value of each cell (NaN outside the boxes)
struct Grid
{
    int W, H;
    double dx, dy;
    int nboxes, nx[SNAP_MAX_BOXES], ny[SNAP_MAX_BOXES], oi[SNAP_MAX_BOXES], oj[SNAP_MAX_BOXES];
    float *val;  // field (log10 if requested)
    float *disp; // displayed value: val, or schlieren intensity
};

static int grid_init(struct Grid *g, const struct Snapshot *s)
{
    g->W = g->H = 0;
    g->dx = s->dx;
    g->dy = s->dy;
    g->nboxes = s->nboxes;
    for (int b = 0; b < s->nboxes; ++b)
    {
        g->nx[b] = s->nx[b];
        g->ny[b] = s->ny[b];
        g->oi[b] = s->origin_i[b];
        g->oj[b] = s->origin_j[b];
        if (g->oi[b] < 0 || g->oj[b] < 0)
            return -1;
        if (g->oi[b] + g->nx[b] > g->W)
            g->W = g->oi[b] + g->nx[b];
        if (g->oj[b] + g->ny[b] > g->H)
            g->H = g->oj[b] + g->ny[b];
    }
    size_t n = (size_t)g->W * g->H;
    g->val = malloc(n * sizeof(float));
    g->disp = malloc(n * sizeof(float));
    if (!g->val || !g->disp)
        return -1;
    for (size_t k = 0; k < n; ++k)
        g->val[k] = g->disp[k] = NAN;
    return 0;
}

static int grid_matches(const struct Grid *g, const struct Snapshot *s)
{
    if (s->nboxes != g->nboxes)
        return 0;
    for (int b = 0; b < s->nboxes; ++b)
        if (s->nx[b] != g->nx[b] || s->ny[b] != g->ny[b] || s->origin_i[b] != g->oi[b] || s->origin_j[b] != g->oj[b])
            return 0;
    return 1;
}

// Fill the grid from a snapshot; returns the maximum schlieren gradient (or 0)
static double grid_fill(struct Grid *g, const struct Snapshot *s, int slot, const struct Options *opt)
{
    for (int b = 0; b < s->nboxes; ++b)
    {
        const float *src = snapshot_box_field(s, b, slot);
        for (int j = 0; j < s->ny[b]; ++j)
        {
            float *dst = g->val + (size_t)(g->oj[b] + j) * g->W + g->oi[b];
            const float *row = src + (size_t)j * s->nx[b];
            for (int i = 0; i < s->nx[b]; ++i)
                dst[i] = opt->log_scale ? log10f(row[i]) : row[i];
        }
    }
    if (!opt->schlieren)
    {
        memcpy(g->disp, g->val, (size_t)g->W * g->H * sizeof(float));
        return 0.0;
    }

    // |grad rho| across box seams; one-sided next to walls and the background
    double gmax = 0.0;
    const int W = g->W, H = g->H;
    for (int j = 0; j < H; ++j)
        for (int i = 0; i < W; ++i)
        {
            size_t k = (size_t)j * W + i;
            float c = g->val[k];
            if (!isfinite(c))
                continue;
            float l = i > 0 ? g->val[k - 1] : NAN, r = i < W - 1 ? g->val[k + 1] : NAN;
            float d = j > 0 ? g->val[k - W] : NAN, u = j < H - 1 ? g->val[k + W] : NAN;
            double gx = isfinite(l) && isfinite(r) ? (r - l) / (2.0 * g->dx)
                        : isfinite(r)              ? (r - c) / g->dx
                        : isfinite(l)              ? (c - l) / g->dx
                                                   : 0.0;
            double gy = isfinite(d) && isfinite(u) ? (u - d) / (2.0 * g->dy)
                        : isfinite(u)              ? (u - c) / g->dy
                        : isfinite(d)              ? (c - d) / g->dy
                                                   : 0.0;
            double grad = sqrt(gx * gx + gy * gy);
            g->disp[k] = (float)grad;
            if (grad > gmax)
                gmax = grad;
        }
    return gmax;
}

static void grid_minmax(const struct Grid *g, double *lo, double *hi)
{
    size_t n = (size_t)g->W * g->H;
    for (size_t k = 0; k < n; ++k)
    {
        float v = g->disp[k];
        if (isfinite(v))
        {
            if (v < *lo)
                *lo = v;
            if (v > *hi)
                *hi = v;
        }
    }
}

static int even(double x)
{
    int n = (int)(x + 0.5);
    n += n & 1;
    return n < 2 ? 2 : n;
}

int main(int argc, char **argv)
{
    struct Options opt = {0};
    opt.ffmpeg = "ffmpeg";
    opt.preset = "slow";
    opt.field = SNAP_FIELD_RHO;
    opt.schlieren_k = 15.0;
    opt.scale = 1.0;
    opt.physical = 1;
    opt.fps = 30;
    opt.crf = 26;
    opt.stride = 1;
    opt.still = -1;
    opt.bg[0] = opt.bg[1] = opt.bg[2] = 40;
    opt.gray = -1;
    if (parse_args(argc, argv, &opt) != 0)
    {
        usage(argv[0]);
        return 1;
    }
    // Schlieren images are gray and always computed from the density
    if (opt.gray < 0)
        opt.gray = opt.schlieren;
    if (opt.schlieren)
        opt.field = SNAP_FIELD_RHO;

    // Frames to use: all of them, or just one for --still
    int nfiles = 0;
    char **files = list_frames(opt.dir, &nfiles);
    if (!files || nfiles == 0)
    {
        fprintf(stderr, "No frame_*.snap files in %s\n", opt.dir);
        return 1;
    }
    int first = 0, last = nfiles - 1;
    if (opt.still >= 0)
    {
        if (opt.still >= nfiles)
        {
            fprintf(stderr, "Frame %d out of range (%d frames)\n", opt.still, nfiles);
            return 1;
        }
        first = last = opt.still;
    }
    if (opt.out[0] == '\0')
        snprintf(opt.out, sizeof(opt.out), "%s/%s%s.%s", opt.dir, opt.schlieren ? "schlieren_" : "",
                 snapshot_field_name(opt.field), opt.still >= 0 ? "ppm" : "mp4");

    // The first frame defines the geometry; every other frame must match it
    struct Snapshot s = {0};
    struct Grid g = {0};
    if (snapshot_read(files[first], &s) != 0 || grid_init(&g, &s) != 0)
        return 1;
    int slot = snapshot_field_slot(&s, opt.field);
    if (slot < 0)
    {
        fprintf(stderr, "Field %s is not in the snapshots (write it with --fields=...)\n",
                snapshot_field_name(opt.field));
        return 1;
    }

    // Pass 1: global range over the selected frames
    double lo = INFINITY, hi = -INFINITY, gmax = 0.0;
    int need_pass1 = opt.schlieren || !opt.have_range;
    for (int f = first; f <= last && need_pass1; f += opt.stride)
    {
        if (snapshot_read(files[f], &s) != 0 || !grid_matches(&g, &s) || snapshot_field_slot(&s, opt.field) != slot)
        {
            fprintf(stderr, "Inconsistent frame %s\n", files[f]);
            return 1;
        }
        double m = grid_fill(&g, &s, slot, &opt);
        if (m > gmax)
            gmax = m;
        if (!opt.schlieren)
            grid_minmax(&g, &lo, &hi);
    }
    if (opt.schlieren)
    {
        lo = 0.0;
        hi = 1.0;
    }
    else if (opt.have_range)
    {
        lo = opt.lo;
        hi = opt.hi;
    }
    if (!(hi > lo))
        hi = lo + 1.0;
    const double inv = 1.0 / (hi - lo);

    // Output size and pixel -> cell map (row 0 of the image is the top of the domain)
    double aspect = opt.physical ? (g.W * g.dx) / (g.H * g.dy) : (double)g.W / g.H;
    const int Hout = even(g.H * opt.scale);
    const int Wout = even(Hout * aspect);
    size_t npix = (size_t)Wout * Hout;
    size_t *pixmap = malloc(npix * sizeof(size_t));
    unsigned char *rgb = malloc(npix * 3);
    if (!pixmap || !rgb)
        return 1;
    for (int py = 0; py < Hout; ++py)
    {
        int gj = g.H - 1 - (int)((py + 0.5) * g.H / Hout);
        for (int px = 0; px < Wout; ++px)
        {
            int gi = (int)((px + 0.5) * g.W / Wout);
            pixmap[(size_t)py * Wout + px] = (size_t)gj * g.W + gi;
        }
    }
    unsigned char lut[256][3];
    build_lut(lut, opt.gray);

    FILE *out;
    if (opt.still >= 0)
        out = fopen(opt.out, "wb");
    else
    {
        char cmd[2048];
        snprintf(cmd, sizeof(cmd),
                 "\"%s\" -loglevel error -y -f rawvideo -pix_fmt rgb24 -s %dx%d -r %d -i - "
                 "-c:v libx264 -preset %s -crf %d -pix_fmt yuv420p -movflags +faststart \"%s\"",
                 opt.ffmpeg, Wout, Hout, opt.fps, opt.preset, opt.crf, opt.out);
        out = popen(cmd, "w");
    }
    if (!out)
    {
        perror(opt.still >= 0 ? opt.out : opt.ffmpeg);
        return 1;
    }
    if (opt.still >= 0)
        fprintf(out, "P6\n%d %d\n255\n", Wout, Hout);

    fprintf(stderr, "%s: %d frame(s), %dx%d, %s%s range [%g, %g] -> %s\n", opt.dir, (last - first) / opt.stride + 1,
            Wout, Hout, opt.schlieren ? "schlieren " : "", snapshot_field_name(opt.field), lo, hi, opt.out);

    // Pass 2: compose and encode
    int status = 0, written = 0;
    for (int f = first; f <= last; f += opt.stride)
    {
        if (snapshot_read(files[f], &s) != 0 || !grid_matches(&g, &s))
        {
            status = 1;
            break;
        }
        grid_fill(&g, &s, slot, &opt);
        for (size_t p = 0; p < npix; ++p)
        {
            double v = g.disp[pixmap[p]];
            const unsigned char *c = opt.bg;
            if (isfinite(v))
            {
                if (opt.schlieren)
                    v = exp(-opt.schlieren_k * v / (gmax > 0.0 ? gmax : 1.0));
                double t = (v - lo) * inv;
                t = t < 0.0 ? 0.0 : t > 1.0 ? 1.0 : t;
                c = lut[(int)(255.0 * t + 0.5)];
            }
            rgb[3 * p] = c[0];
            rgb[3 * p + 1] = c[1];
            rgb[3 * p + 2] = c[2];
        }
        if (fwrite(rgb, 1, npix * 3, out) != npix * 3)
        {
            fprintf(stderr, "Error writing frame %d\n", f);
            status = 1;
            break;
        }
        if (++written % 100 == 0)
            fprintf(stderr, "  %d frames\n", written);
    }

    int rc = opt.still >= 0 ? fclose(out) : pclose(out);
    if (rc != 0)
    {
        fprintf(stderr, "Error: %s failed (is ffmpeg installed? try --ffmpeg=PATH)\n", opt.still >= 0 ? "write" : "ffmpeg");
        status = 1;
    }

    for (int k = 0; k < nfiles; ++k)
        free(files[k]);
    free(files);
    free(pixmap);
    free(rgb);
    free(g.val);
    free(g.disp);
    snapshot_free(&s);
    return status;
}
