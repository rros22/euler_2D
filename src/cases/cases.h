// Problem setups ("cases"). A case builds the boxes of the domain, sets their boundary conditions
// and initial data, and records where each box sits in a global cell grid (used for output and
// rendering only). To add a case: write a setup function like the ones below (see step.c for a
// short example), declare it here and add it to the table in cases.c.
#ifndef CASES_H
#define CASES_H

#include <stdio.h>
#include "mesh/box.h"

#define CASE_MAX_BOXES 8

struct Case
{
    const char *name;
    int num_boxes;
    struct Box *boxes[CASE_MAX_BOXES];
    const char *box_names[CASE_MAX_BOXES];
    int origin_i[CASE_MAX_BOXES], origin_j[CASE_MAX_BOXES]; // global cell of interior cell (1, 1), 0-based
    double dx, dy, gamma;                                   // one cell size for all boxes
    double t_end;                                           // default end time
    int max_steps;                                          // default step limit
};

// Helpers for writing cases

// Start a case: clears c and sets the name, cell size, gamma and default run length
void case_begin(struct Case *c, const char *name, double dx, double dy, double gamma, double t_end, int max_steps);

// Add a box of nx x ny cells whose lower-left cell sits at global cell (origin_i, origin_j).
// Returns NULL after printing an error if allocation fails or the case already has CASE_MAX_BOXES.
struct Box *case_add_box(struct Case *c, const char *name, int nx, int ny, int origin_i, int origin_j);

// Glue edge_a of box a to edge_b of box b, in both directions (e.g. EDGE_TOP of a to EDGE_BOTTOM of b)
void case_connect(struct Box *a, int edge_a, struct Box *b, int edge_b);

// Last step of every case: compute the conserved variables from the primitive initial data
void case_finish(struct Case *c);

void case_destroy(struct Case *c);

// The cases. base sets the resolution; its meaning is case specific (see each file).

// Original 2D jet: a pressurized chamber (top box) discharging through the centre bottom box into
// an ambient region (left, bottom, right boxes). base = resolution unit (original: 100).
int case_jet(struct Case *c, int base, double gamma);

// Sod shock tube with n cells along x (dir = 0) or y (dir = 1) and 4 cells across
int case_sod(struct Case *c, int n, int dir, double gamma);

// Mach 3 wind tunnel with a forward-facing step; base = cells per 0.2 length (default 20)
int case_step(struct Case *c, int base, double gamma);

// Dense subsonic jet injected from the top into lighter gas at rest; base = cells per nozzle width (default 48)
int case_injection(struct Case *c, int base, double gamma);

// Build a case by name with base cells (base <= 0: the case's default); returns 0 on success
int case_create(struct Case *c, const char *name, int base, double gamma);

// Print the available cases, one per line, for --help
void case_print_list(FILE *fp);

#endif
