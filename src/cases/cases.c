// Case registry and the helpers shared by all case setups.
// Every case in CASES[] needs a section in docs/cases.md (checked by make test).
#include "cases/cases.h"
#include <string.h>

// Sod along x, so it fits the common (case, base, gamma) signature
static int case_sod_x(struct Case *c, int n, double gamma)
{
    return case_sod(c, n, 0, gamma);
}

struct CaseInfo
{
    const char *name;
    int (*setup)(struct Case *c, int base, double gamma);
    int default_base;
    const char *description;
};

// To add a case, add a line here
static const struct CaseInfo CASES[] = {
    {"jet", case_jet, 100, "original 2D jet: pressurized chamber discharging into ambient gas (base 100: 500x800)"},
    {"sod", case_sod_x, 200, "Sod shock tube along x (base = number of cells)"},
    {"step", case_step, 20, "Mach 3 wind tunnel with a forward-facing step (base 20: 300x100)"},
    {"injection", case_injection, 48, "dense subsonic jet injected into lighter gas (base 48: 432x576)"},
};
static const int NUM_CASES = sizeof(CASES) / sizeof(CASES[0]);

int case_create(struct Case *c, const char *name, int base, double gamma)
{
    for (int k = 0; k < NUM_CASES; ++k)
        if (strcmp(name, CASES[k].name) == 0)
            return CASES[k].setup(c, base > 0 ? base : CASES[k].default_base, gamma);

    fprintf(stderr, "Unknown case '%s'. Available cases:\n", name);
    case_print_list(stderr);
    return -1;
}

void case_print_list(FILE *fp)
{
    for (int k = 0; k < NUM_CASES; ++k)
        fprintf(fp, "    %-10s %s\n", CASES[k].name, CASES[k].description);
}

void case_begin(struct Case *c, const char *name, double dx, double dy, double gamma, double t_end, int max_steps)
{
    memset(c, 0, sizeof(*c));
    c->name = name;
    c->dx = dx;
    c->dy = dy;
    c->gamma = gamma;
    c->t_end = t_end;
    c->max_steps = max_steps;
}

struct Box *case_add_box(struct Case *c, const char *name, int nx, int ny, int origin_i, int origin_j)
{
    if (c->num_boxes == CASE_MAX_BOXES)
    {
        fprintf(stderr, "Error: case %s has more than %d boxes\n", c->name, CASE_MAX_BOXES);
        return NULL;
    }
    struct Box *b = box_create(nx, ny);
    if (!b)
    {
        fprintf(stderr, "Error: could not allocate box %s (%d x %d) of case %s\n", name, nx, ny, c->name);
        return NULL;
    }
    int k = c->num_boxes++;
    c->boxes[k] = b;
    c->box_names[k] = name;
    c->origin_i[k] = origin_i;
    c->origin_j[k] = origin_j;
    return b;
}

void case_connect(struct Box *a, int edge_a, struct Box *b, int edge_b)
{
    box_set_boundary_connected(a, edge_a, b, edge_b);
    box_set_boundary_connected(b, edge_b, a, edge_a);
}

void case_finish(struct Case *c)
{
    for (int k = 0; k < c->num_boxes; ++k)
        box_prim_to_cons(c->boxes[k], c->gamma);
}

void case_destroy(struct Case *c)
{
    for (int k = 0; k < c->num_boxes; ++k)
        box_destroy(c->boxes[k]);
    memset(c, 0, sizeof(*c));
}
