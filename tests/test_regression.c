// Bitwise regression of the solver against stored golden states (first order and MUSCL-Hancock).
// Usage: test_regression [--write]
#include "test_util.h"
#include "cases/cases.h"
#include <stdio.h>
#include <string.h>

#define GAMMA 1.4
#define CFL 0.4
#define BASE 10
#define STEPS 40

static int run_case(const char *path, const struct SchemeConfig *cfg, int write)
{
    struct Setup s;
    if (setup_four_box(&s, BASE, 0, GAMMA) != 0)
    {
        fprintf(stderr, "setup failed\n");
        return 1;
    }
    setup_run(&s, GAMMA, CFL, 1e30, STEPS, cfg);

    int status = 0;
    if (write)
    {
        status = setup_dump(&s, path) == 0 ? 0 : 1;
        printf("%s golden %s\n", status ? "FAILED to write" : "wrote", path);
    }
    else
    {
        double max_diff;
        long ndiff = setup_compare(&s, path, &max_diff);
        if (ndiff < 0)
        {
            printf("FAIL regression: cannot read %s\n", path);
            status = 1;
        }
        else if (ndiff > 0)
        {
            printf("FAIL regression %s: %ld values differ, max |diff| = %.3e\n", path, ndiff, max_diff);
            status = 1;
        }
        else
            printf("PASS regression: bitwise identical to %s\n", path);
    }

    setup_destroy(&s);
    return status;
}

// The jet case must reproduce the setup the goldens were made from (the original main.c geometry)
static int run_jet_case(const char *path, const struct SchemeConfig *cfg)
{
    struct Case c;
    if (case_jet(&c, BASE, GAMMA) != 0 || c.num_boxes > SETUP_MAX_BOXES)
    {
        fprintf(stderr, "case_jet failed\n");
        return 1;
    }
    struct Setup s = {0};
    s.n = c.num_boxes;
    s.dx = c.dx;
    s.dy = c.dy;
    for (int b = 0; b < c.num_boxes; ++b)
        s.boxes[b] = c.boxes[b];
    setup_run(&s, GAMMA, CFL, 1e30, STEPS, cfg);

    double max_diff;
    long ndiff = setup_compare(&s, path, &max_diff);
    int status = ndiff != 0;
    if (status)
        printf("FAIL jet case vs %s: %ld values differ, max |diff| = %.3e\n", path, ndiff, max_diff);
    else
        printf("PASS jet case: bitwise identical to %s\n", path);
    case_destroy(&c);
    return status;
}

int main(int argc, char **argv)
{
    int write = argc > 1 && strcmp(argv[1], "--write") == 0;
    const struct SchemeConfig first = {RECON_FIRST_ORDER, LIMITER_MINMOD};
    const struct SchemeConfig muscl = {RECON_MUSCL_HANCOCK, LIMITER_VAN_LEER};

    int status = run_case("tests/golden/first_order_4box.bin", &first, write);
    status |= run_case("tests/golden/muscl_vanleer_4box.bin", &muscl, write);
    if (!write)
    {
        status |= run_jet_case("tests/golden/first_order_4box.bin", &first);
        status |= run_jet_case("tests/golden/muscl_vanleer_4box.bin", &muscl);
    }
    return status;
}
