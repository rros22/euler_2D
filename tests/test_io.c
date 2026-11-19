// Snapshot round trip: write the four-box setup after a few steps, read it back, compare with the boxes.
// Usage: test_io [scratch_dir]   (default: build)
#include "test_util.h"
#include "io/snapshot.h"
#include <stdio.h>
#include <string.h>

#define GAMMA 1.4

static int failures = 0;

static void check(int ok, const char *what)
{
    if (!ok)
    {
        printf("FAIL io: %s\n", what);
        ++failures;
    }
}

int main(int argc, char **argv)
{
    char path[512];
    snprintf(path, sizeof(path), "%s/test_io.snap", argc > 1 ? argv[1] : "build");

    // Field list parsing
    int f[SNAP_MAX_FIELDS];
    check(snapshot_parse_fields("rho", f) == 1 && f[0] == SNAP_FIELD_RHO, "parse rho");
    check(snapshot_parse_fields("p,u,v,rho", f) == 4 && f[0] == SNAP_FIELD_P && f[3] == SNAP_FIELD_RHO, "parse list");
    check(snapshot_parse_fields("rho,rho", f) < 0, "reject duplicate field");
    check(snapshot_parse_fields("rho,e", f) < 0, "reject unknown field");
    check(snapshot_parse_fields("", f) < 0, "reject empty list");

    struct Setup s;
    if (setup_four_box(&s, 6, 0, GAMMA) != 0)
        return 1;
    const struct SchemeConfig cfg = {RECON_MUSCL_HANCOCK, LIMITER_VAN_LEER};
    double t = setup_run(&s, GAMMA, 0.4, 1e30, 5, &cfg);

    int oi[4] = {0, 12, 18, 12}, oj[4] = {0, 0, 0, 36};
    int fields[3] = {SNAP_FIELD_P, SNAP_FIELD_RHO, SNAP_FIELD_V};
    check(snapshot_write(path, 5, t, s.dx, s.dy, s.boxes, s.n, oi, oj, fields, 3) == 0, "write");

    struct Snapshot r = {0};
    if (snapshot_read(path, &r) != 0)
    {
        printf("FAIL io: read %s\n", path);
        return 1;
    }
    check(r.step == 5 && r.time == t && r.dx == s.dx && r.dy == s.dy, "header step/time/dx/dy");
    check(r.nfields == 3 && r.fields[0] == SNAP_FIELD_P && r.fields[1] == SNAP_FIELD_RHO && r.fields[2] == SNAP_FIELD_V,
          "header fields");
    check(r.nboxes == s.n, "box count");
    check(snapshot_field_slot(&r, SNAP_FIELD_RHO) == 1 && snapshot_field_slot(&r, SNAP_FIELD_U) == -1, "field slots");

    long mismatches = 0;
    for (int b = 0; b < r.nboxes && b < s.n; ++b)
    {
        const struct Box *box = s.boxes[b];
        check(r.nx[b] == box->nx && r.ny[b] == box->ny && r.origin_i[b] == oi[b] && r.origin_j[b] == oj[b],
              "box geometry");
        const double *src[3] = {box->w_4, box->w_1, box->w_3};
        for (int k = 0; k < 3; ++k)
        {
            const float *v = snapshot_box_field(&r, b, k);
            for (int j = 1; j <= box->ny; ++j)
                for (int i = 1; i <= box->nx; ++i)
                    mismatches += v[(j - 1) * box->nx + (i - 1)] != (float)src[k][INDEX(i, j, box->nx)];
        }
    }
    check(mismatches == 0, "values equal to float of the box data");

    // Corrupt magic must be rejected (snapshot_read reports the error on stderr)
    fprintf(stderr, "(expected error follows)\n");
    FILE *fp = fopen(path, "r+b");
    if (fp)
    {
        fputc('X', fp);
        fclose(fp);
    }
    check(snapshot_read(path, &r) != 0, "reject bad magic");
    remove(path);

    snapshot_free(&r);
    setup_destroy(&s);
    printf("%s io: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
    return failures != 0;
}
