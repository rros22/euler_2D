#ifndef TEST_UTIL_H
#define TEST_UTIL_H

#include "mesh/box.h"
#include "numerics/sweep.h"

#define SETUP_MAX_BOXES 4

struct Setup
{
    int n;
    struct Box *boxes[SETUP_MAX_BOXES];
    double dx, dy;
};

// Geometry and initial data of main.c, scaled by base; closed = 1 makes every outer edge reflective
int setup_four_box(struct Setup *s, int base, int closed, double gamma);

// Single box with reflective edges; caller sets boundaries and primitive initial data, then calls setup_finish
int setup_single_box(struct Setup *s, int nx, int ny, double dx, double dy);
void setup_finish(struct Setup *s, double gamma);
void setup_destroy(struct Setup *s);

// Advance until t_end (last step shortened to land on it) or max_steps; returns the reached time
double setup_run(struct Setup *s, double gamma, double cfl, double t_end, int max_steps,
                 const struct SchemeConfig *cfg);

// Sum of a conserved variable (k = 1..4) over interior cells of all boxes, times cell area
double setup_total(const struct Setup *s, int k);

// 1 if all interior rho and p are finite and positive
int setup_is_physical(const struct Setup *s);

// Full-precision dump of interior w and u of all boxes
int setup_dump(const struct Setup *s, const char *path);

// Returns number of differing doubles (-1 on error); max_abs_diff receives the largest difference
long setup_compare(const struct Setup *s, const char *path, double *max_abs_diff);

#endif
