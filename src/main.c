// euler_2D: 2D compressible Euler solver (finite volume, exact Riemann solver, dimensional splitting)
// on a domain made of connected rectangular boxes.
//
// Run "euler_2D --help" for the options. The program builds a case (src/cases), advances it to the
// end time and writes binary frames at regular simulation times for tools/render.
// Options are documented in docs/options.md; keep usage() and the docs in sync (checked by make test).
#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/stat.h>
#include "cases/cases.h"
#include "io/snapshot.h"
#include "mesh/box.h"
#include "solver/stepper.h"

#define GAMMA 1.4 // ratio of specific heats (air)

// Command line options. Negative values mean "use the case's default".
struct Options
{
    const char *case_name;
    int base;      // resolution, meaning depends on the case
    double t_end;  // end time
    int max_steps; // step limit
    double out_dt; // time between frames; 0 = no frames
    char out_dir[512];
    int fields[SNAP_MAX_FIELDS], nfields; // variables stored in each frame
    int csv;                              // also write CSV files at the start and the end
    double cfl;
    struct SchemeConfig scheme;
};

static void usage(FILE *fp, const char *prog)
{
    fprintf(fp,
            "Usage: %s [options]\n"
            "  --case=NAME               problem to run (default jet), one of:\n",
            prog);
    case_print_list(fp);
    fprintf(fp,
            "  --base=N                  resolution (see the case list; default per case)\n"
            "  --t-end=T                 end time (default per case)\n"
            "  --max-steps=N             step limit (default per case)\n"
            "  --out-dt=DT               time between frames (default t_end/200, 0 = no frames)\n"
            "  --out-dir=DIR             output directory (default results/<case>)\n"
            "  --fields=rho,u,v,p        variables stored in frames (default rho)\n"
            "  --csv                     also write initial/final CSV files per box\n"
            "  --scheme=first|muscl      first order or second order MUSCL-Hancock (default first)\n"
            "  --limiter=minmod|vanleer|mc  slope limiter for MUSCL (default vanleer)\n"
            "  --cfl=C                   CFL number (default 0.4)\n"
            "  --help                    show this message\n");
}

static void set_defaults(struct Options *opt)
{
    memset(opt, 0, sizeof(*opt));
    opt->case_name = "jet";
    opt->t_end = -1.0;
    opt->max_steps = -1;
    opt->out_dt = -1.0;
    opt->fields[0] = SNAP_FIELD_RHO;
    opt->nfields = 1;
    opt->cfl = 0.4;
    opt->scheme.recon = RECON_FIRST_ORDER;
    opt->scheme.limiter = LIMITER_VAN_LEER;
}

// Returns 0 to run, 1 if the program should exit successfully (--help), -1 on a bad option
static int parse_args(int argc, char **argv, struct Options *opt)
{
    for (int a = 1; a < argc; ++a)
    {
        const char *arg = argv[a];
        const char *val = strchr(arg, '=');
        val = val ? val + 1 : "";

        if (strcmp(arg, "--help") == 0 || strcmp(arg, "-h") == 0)
        {
            usage(stdout, argv[0]);
            return 1;
        }
        else if (strncmp(arg, "--case=", 7) == 0)
            opt->case_name = val;
        else if (strncmp(arg, "--base=", 7) == 0)
            opt->base = atoi(val);
        else if (strncmp(arg, "--t-end=", 8) == 0)
            opt->t_end = atof(val);
        else if (strncmp(arg, "--max-steps=", 12) == 0)
            opt->max_steps = atoi(val);
        else if (strncmp(arg, "--out-dt=", 9) == 0)
            opt->out_dt = atof(val);
        else if (strncmp(arg, "--out-dir=", 10) == 0)
            snprintf(opt->out_dir, sizeof(opt->out_dir), "%s", val);
        else if (strncmp(arg, "--fields=", 9) == 0)
        {
            opt->nfields = snapshot_parse_fields(val, opt->fields);
            if (opt->nfields < 0)
            {
                fprintf(stderr, "Invalid --fields '%s' (choose from rho,u,v,p)\n", val);
                return -1;
            }
        }
        else if (strcmp(arg, "--csv") == 0)
            opt->csv = 1;
        else if (strcmp(arg, "--scheme=first") == 0)
            opt->scheme.recon = RECON_FIRST_ORDER;
        else if (strcmp(arg, "--scheme=muscl") == 0)
            opt->scheme.recon = RECON_MUSCL_HANCOCK;
        else if (strcmp(arg, "--limiter=minmod") == 0)
            opt->scheme.limiter = LIMITER_MINMOD;
        else if (strcmp(arg, "--limiter=vanleer") == 0)
            opt->scheme.limiter = LIMITER_VAN_LEER;
        else if (strcmp(arg, "--limiter=mc") == 0)
            opt->scheme.limiter = LIMITER_MC;
        else if (strncmp(arg, "--cfl=", 6) == 0)
            opt->cfl = atof(val);
        else
        {
            fprintf(stderr, "Unknown option %s\n", arg);
            usage(stderr, argv[0]);
            return -1;
        }
    }
    return 0;
}

// Create a directory and its parents (like mkdir -p)
static int make_dirs(const char *path)
{
    char buf[512];
    snprintf(buf, sizeof(buf), "%s", path);
    for (char *p = buf + 1;; ++p)
    {
        if (*p == '/' || *p == '\0')
        {
            char c = *p;
            *p = '\0';
            if (mkdir(buf, 0755) != 0 && errno != EEXIST)
            {
                perror(buf);
                return -1;
            }
            *p = c;
            if (c == '\0')
                return 0;
        }
    }
}

static const char *limiter_name(enum Limiter l)
{
    return l == LIMITER_MINMOD ? "minmod" : l == LIMITER_MC ? "mc" : "vanleer";
}

// case.txt: a human-readable record of the run settings and the box layout (not read by any tool)
static void write_case_info(const struct Options *opt, const struct Case *problem)
{
    char path[600];
    snprintf(path, sizeof(path), "%s/case.txt", opt->out_dir);
    FILE *fp = fopen(path, "w");
    if (!fp)
    {
        perror(path);
        return;
    }
    fprintf(fp, "case %s\nscheme %s\nlimiter %s\ncfl %g\ngamma %g\ndx %.17g\ndy %.17g\nt_end %g\nout_dt %g\n",
            problem->name, opt->scheme.recon == RECON_MUSCL_HANCOCK ? "muscl" : "first",
            limiter_name(opt->scheme.limiter), opt->cfl, problem->gamma, problem->dx, problem->dy, opt->t_end,
            opt->out_dt);
    fprintf(fp, "fields");
    for (int f = 0; f < opt->nfields; ++f)
        fprintf(fp, " %s", snapshot_field_name(opt->fields[f]));
    fprintf(fp, "\nboxes %d\n", problem->num_boxes);
    for (int b = 0; b < problem->num_boxes; ++b)
        fprintf(fp, "box %s nx %d ny %d origin_i %d origin_j %d\n", problem->box_names[b], problem->boxes[b]->nx,
                problem->boxes[b]->ny, problem->origin_i[b], problem->origin_j[b]);
    fclose(fp);
}

// One CSV per box (box_<name>_<tag>.csv) plus the connectivity table
static void write_csv(const struct Options *opt, struct Case *problem, const char *tag)
{
    char path[700];
    for (int b = 0; b < problem->num_boxes; ++b)
    {
        snprintf(path, sizeof(path), "%s/box_%s_%s.csv", opt->out_dir, problem->box_names[b], tag);
        save_box_to_csv(problem->boxes[b], path);
    }
    snprintf(path, sizeof(path), "%s/box_connections.csv", opt->out_dir);
    save_box_connections(problem->boxes, problem->num_boxes, path);
}

static int write_frame(const struct Options *opt, const struct Case *problem, int frame, int step, double t)
{
    char path[600];
    snprintf(path, sizeof(path), "%s/frame_%05d.snap", opt->out_dir, frame);
    return snapshot_write(path, step, t, problem->dx, problem->dy, problem->boxes, problem->num_boxes,
                          problem->origin_i, problem->origin_j, opt->fields, opt->nfields);
}

int main(int argc, char **argv)
{
    // Options, then the case, then the defaults that depend on the case
    struct Options opt;
    set_defaults(&opt);
    int parsed = parse_args(argc, argv, &opt);
    if (parsed != 0)
        return parsed > 0 ? 0 : 1;

    struct Case problem;
    if (case_create(&problem, opt.case_name, opt.base, GAMMA) != 0)
        return 1;
    if (opt.t_end < 0.0)
        opt.t_end = problem.t_end;
    if (opt.max_steps < 0)
        opt.max_steps = problem.max_steps;
    if (opt.out_dt < 0.0)
        opt.out_dt = opt.t_end / 200.0;
    if (opt.out_dir[0] == '\0')
        snprintf(opt.out_dir, sizeof(opt.out_dir), "results/%s", problem.name);

    if (make_dirs(opt.out_dir) != 0)
    {
        case_destroy(&problem);
        return 1;
    }
    write_case_info(&opt, &problem);
    if (opt.csv)
        write_csv(&opt, &problem, "initial");

    printf("case %s, scheme %s%s, cfl %g, t_end %g, frames every %g -> %s\n", problem.name,
           opt.scheme.recon == RECON_MUSCL_HANCOCK ? "muscl/" : "first",
           opt.scheme.recon == RECON_MUSCL_HANCOCK ? limiter_name(opt.scheme.limiter) : "", opt.cfl, opt.t_end,
           opt.out_dt, opt.out_dir);

    // Time loop. Frames are written at t = 0, out_dt, 2 out_dt, ...; the time step is shortened when
    // needed so the solution lands exactly on each frame time and on t_end.
    double t = 0.0;
    int step = 0, frame = 0, status = 0;
    double t_next_frame = opt.out_dt > 0.0 ? 0.0 : INFINITY;
    double vmax = solver_max_wave_speed(problem.boxes, problem.num_boxes, GAMMA);
    clock_t wall_start = clock();

    for (;;)
    {
        if (t >= t_next_frame * (1.0 - 1e-12))
        {
            if (write_frame(&opt, &problem, frame, step, t) != 0)
            {
                status = 1;
                break;
            }
            printf("frame %5d  step %6d  t = %.5f\n", frame, step, t);
            ++frame;
            t_next_frame = frame * opt.out_dt;
        }
        if (step >= opt.max_steps || t >= opt.t_end * (1.0 - 1e-12))
            break;
        if (!(vmax > 0.0) || !isfinite(vmax))
        {
            fprintf(stderr, "Error: invalid wave speed %g at step %d (t = %g)\n", vmax, step, t);
            status = 1;
            break;
        }

        // CFL time step. One X-Y-Y-X step advances the solution by 2 * dt.
        double dt = opt.cfl * fmin(problem.dx, problem.dy) / vmax;
        double t_stop = fmin(opt.t_end, t_next_frame);
        if (t + 2.0 * dt > t_stop)
            dt = 0.5 * (t_stop - t);

        vmax = solver_step(problem.boxes, problem.num_boxes, dt, problem.dx, problem.dy, GAMMA, &opt.scheme);
        t += 2.0 * dt;
        step++;
    }

    // Final state as a last frame when t_end is not a multiple of out_dt (or max_steps stopped the run)
    if (status == 0 && opt.out_dt > 0.0 && t > (frame - 1) * opt.out_dt * (1.0 + 1e-12))
    {
        if (write_frame(&opt, &problem, frame, step, t) == 0)
            ++frame;
        else
            status = 1;
    }

    double wall_seconds = (double)(clock() - wall_start) / CLOCKS_PER_SEC;
    printf("done: %d steps, t = %.5f, %d frames, %.1f s\n", step, t, frame, wall_seconds);
    if (opt.csv)
        write_csv(&opt, &problem, "final");

    case_destroy(&problem);
    return status;
}
