// Simulation output: compact binary snapshots (one file per frame, read by tools/render.c and
// tools/snapshot.py) and the original per-box CSV writer. The file layout is described in snapshot.h.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "io/snapshot.h"

static const char SNAP_MAGIC[8] = {'E', 'U', '2', 'D', 'S', 'N', 'A', 'P'};
static const char *const FIELD_NAMES[SNAP_MAX_FIELDS] = {"rho", "u", "v", "p"}; // indexed by enum SnapField

const char *snapshot_field_name(int field)
{
    return (field >= 0 && field < SNAP_MAX_FIELDS) ? FIELD_NAMES[field] : "?";
}

int snapshot_parse_fields(const char *list, int *fields)
{
    int n = 0;
    const char *p = list;
    while (*p)
    {
        size_t len = strcspn(p, ",");
        int id = -1;
        for (int f = 0; f < SNAP_MAX_FIELDS; ++f)
            if (strlen(FIELD_NAMES[f]) == len && strncmp(p, FIELD_NAMES[f], len) == 0)
                id = f;
        if (id < 0 || n == SNAP_MAX_FIELDS)
            return -1;
        for (int k = 0; k < n; ++k)
            if (fields[k] == id)
                return -1;
        fields[n++] = id;
        p += len;
        if (*p == ',')
            ++p;
    }
    return n > 0 ? n : -1;
}

// Integers are stored as exactly 32 bits whatever the size of int
static int write_i32(FILE *fp, int v)
{
    int32_t x = (int32_t)v;
    return fwrite(&x, sizeof(x), 1, fp) == 1 ? 0 : -1;
}

static int read_i32(FILE *fp, int *v)
{
    int32_t x;
    if (fread(&x, sizeof(x), 1, fp) != 1)
        return -1;
    *v = (int)x;
    return 0;
}

int snapshot_write(const char *path, int step, double time, double dx, double dy, struct Box *const *boxes,
                   int nboxes, const int *origin_i, const int *origin_j, const int *fields, int nfields)
{
    if (nboxes < 1 || nboxes > SNAP_MAX_BOXES || nfields < 1 || nfields > SNAP_MAX_FIELDS)
        return -1;

    // Write to a temporary name and rename, so readers never see a partial frame
    char tmp[1024];
    if (snprintf(tmp, sizeof(tmp), "%s.part", path) >= (int)sizeof(tmp))
        return -1;
    FILE *fp = fopen(tmp, "wb");
    if (!fp)
    {
        perror(tmp);
        return -1;
    }

    int err = 0;
    err |= fwrite(SNAP_MAGIC, 1, sizeof(SNAP_MAGIC), fp) != sizeof(SNAP_MAGIC);
    err |= write_i32(fp, SNAP_VERSION);
    err |= write_i32(fp, step);
    err |= fwrite(&time, sizeof(time), 1, fp) != 1;
    err |= fwrite(&dx, sizeof(dx), 1, fp) != 1;
    err |= fwrite(&dy, sizeof(dy), 1, fp) != 1;
    err |= write_i32(fp, nfields);
    for (int f = 0; f < nfields; ++f)
        err |= write_i32(fp, fields[f]);
    err |= write_i32(fp, nboxes);
    int max_nx = 0;
    for (int b = 0; b < nboxes; ++b)
    {
        err |= write_i32(fp, boxes[b]->nx);
        err |= write_i32(fp, boxes[b]->ny);
        err |= write_i32(fp, origin_i[b]);
        err |= write_i32(fp, origin_j[b]);
        if (boxes[b]->nx > max_nx)
            max_nx = boxes[b]->nx;
    }

    // Data: interior cells only, converted to float row by row (half the size of double, plenty for plots)
    float *row = malloc((size_t)max_nx * sizeof(float));
    if (!row)
        err = 1;
    for (int b = 0; b < nboxes && !err; ++b)
    {
        const struct Box *box = boxes[b];
        for (int f = 0; f < nfields && !err; ++f)
        {
            const double *src = fields[f] == SNAP_FIELD_RHO ? box->w_1
                                : fields[f] == SNAP_FIELD_U ? box->w_2
                                : fields[f] == SNAP_FIELD_V ? box->w_3
                                                            : box->w_4;
            for (int j = 1; j <= box->ny && !err; ++j)
            {
                const double *line = src + INDEX(1, j, box->nx);
                for (int i = 0; i < box->nx; ++i)
                    row[i] = (float)line[i];
                err |= fwrite(row, sizeof(float), (size_t)box->nx, fp) != (size_t)box->nx;
            }
        }
    }
    free(row);

    if (fclose(fp) != 0)
        err = 1;
    if (err || rename(tmp, path) != 0)
    {
        fprintf(stderr, "Error writing snapshot %s\n", path);
        remove(tmp);
        return -1;
    }
    return 0;
}

int snapshot_read(const char *path, struct Snapshot *s)
{
    FILE *fp = fopen(path, "rb");
    if (!fp)
    {
        perror(path);
        return -1;
    }

    // Header, validated as it is read; any failure turns err on and skips the rest
    char magic[8];
    int err = fread(magic, 1, sizeof(magic), fp) != sizeof(magic) || memcmp(magic, SNAP_MAGIC, sizeof(magic)) != 0;
    err = err || read_i32(fp, &s->version) || s->version != SNAP_VERSION;
    err = err || read_i32(fp, &s->step) || fread(&s->time, sizeof(s->time), 1, fp) != 1;
    err = err || fread(&s->dx, sizeof(s->dx), 1, fp) != 1 || fread(&s->dy, sizeof(s->dy), 1, fp) != 1;
    err = err || read_i32(fp, &s->nfields) || s->nfields < 1 || s->nfields > SNAP_MAX_FIELDS;
    for (int f = 0; f < s->nfields && !err; ++f)
        err = read_i32(fp, &s->fields[f]);
    err = err || read_i32(fp, &s->nboxes) || s->nboxes < 1 || s->nboxes > SNAP_MAX_BOXES;

    size_t total = 0;
    for (int b = 0; b < s->nboxes && !err; ++b)
    {
        err = read_i32(fp, &s->nx[b]) || read_i32(fp, &s->ny[b]) || read_i32(fp, &s->origin_i[b]) ||
              read_i32(fp, &s->origin_j[b]) || s->nx[b] < 1 || s->ny[b] < 1;
        s->offset[b] = total;
        total += (size_t)s->nfields * s->nx[b] * s->ny[b];
    }

    // Grow the data buffer only when needed, so reading many frames reuses one allocation
    if (!err && total > s->capacity)
    {
        float *data = realloc(s->data, total * sizeof(float));
        if (data)
        {
            s->data = data;
            s->capacity = total;
        }
        else
            err = 1;
    }
    err = err || fread(s->data, sizeof(float), total, fp) != total;
    fclose(fp);
    if (err)
        fprintf(stderr, "Error: %s is not a valid snapshot\n", path);
    return err ? -1 : 0;
}

void snapshot_free(struct Snapshot *s)
{
    free(s->data);
    s->data = NULL;
    s->capacity = 0;
}

int snapshot_field_slot(const struct Snapshot *s, int field)
{
    for (int f = 0; f < s->nfields; ++f)
        if (s->fields[f] == field)
            return f;
    return -1;
}

// CSV output (main --csv): all interior cells of one box with primitive and conserved variables
void save_box_to_csv(struct Box *box, const char *filename)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("Error opening file");
        return;
    }

    // CSV header
    fprintf(fp, "i,j,w_1,w_2,w_3,w_4,u_1,u_2,u_3,u_4\n");

    // Loop over interior cells only (skip ghost cells)
    for (int j = 1; j <= box->ny; ++j)
    {
        for (int i = 1; i <= box->nx; ++i)
        {
            int idx = INDEX(i, j, box->nx);

            fprintf(fp, "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n",
                    i - 1, j - 1, // shift index to start from 0
                    box->w_1[idx], box->w_2[idx], box->w_3[idx], box->w_4[idx],
                    box->u_1[idx], box->u_2[idx], box->u_3[idx], box->u_4[idx]);
        }
    }

    fclose(fp);
}

// CSV table of every edge's boundary type and, for connected edges, the neighbor box and edge
void save_box_connections(struct Box *boxes[], int count, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (!file)
    {
        fprintf(stderr, "Could not open file %s for writing\n", filename);
        return;
    }

    fprintf(file, "BoxID,Edge,Type,ConnectedToBox,ConnectedToEdge\n");

    for (int b = 0; b < count; ++b)
    {
        struct Box *box = boxes[b];
        for (int edge = 0; edge < 4; ++edge)
        {
            struct BoundaryCondition *bc = &box->boundaries[edge];

            const char *edge_name = (edge == EDGE_LEFT) ? "LEFT" : (edge == EDGE_RIGHT) ? "RIGHT"
                                                               : (edge == EDGE_BOTTOM)  ? "BOTTOM"
                                                                                        : "TOP";

            const char *type_str = (bc->type == BOUNDARY_REFLECTIVE)     ? "REFLECTIVE"
                                   : (bc->type == BOUNDARY_TRANSMISSIVE) ? "TRANSMISSIVE"
                                   : (bc->type == BOUNDARY_INFLOW)       ? "INFLOW"
                                                                         : "CONNECTED";

            int neighbor_id = -1;
            const char *neighbor_edge_name = "";

            if (bc->type == BOUNDARY_CONNECTED)
            {
                for (int n = 0; n < count; ++n)
                {
                    if (boxes[n] == bc->neighbor)
                    {
                        neighbor_id = n;
                        break;
                    }
                }

                neighbor_edge_name = (bc->neighbor_edge == EDGE_LEFT) ? "LEFT" : (bc->neighbor_edge == EDGE_RIGHT) ? "RIGHT"
                                                                             : (bc->neighbor_edge == EDGE_BOTTOM)  ? "BOTTOM"
                                                                                                                   : "TOP";
            }

            fprintf(file, "%d,%s,%s,%d,%s\n", b, edge_name, type_str, neighbor_id, neighbor_edge_name);
        }
    }

    fclose(file);
}
