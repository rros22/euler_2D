#ifndef SNAPSHOT_H
#define SNAPSHOT_H

#include <stddef.h>
#include "mesh/box.h"

// Binary snapshot (native byte order, little-endian on all supported machines):
//   char    magic[8] = "EU2DSNAP"
//   int32   version, step
//   float64 time, dx, dy
//   int32   nfields, field_id[nfields]          (SNAP_FIELD_*)
//   int32   nboxes, then per box: nx, ny, origin_i, origin_j
//   float32 data: for each box, for each field, nx * ny interior values (i fastest, j = 1 at the bottom)
// origin_i/j place interior cell (1, 1) of the box in a global, 0-based cell grid.

#define SNAP_VERSION 1
#define SNAP_MAX_FIELDS 4
#define SNAP_MAX_BOXES 64

enum SnapField
{
    SNAP_FIELD_RHO,
    SNAP_FIELD_U,
    SNAP_FIELD_V,
    SNAP_FIELD_P,
};

const char *snapshot_field_name(int field);

// Parse a comma separated list such as "rho,u,v,p"; returns the number of fields or -1
int snapshot_parse_fields(const char *list, int *fields);

// Write the primitive fields of the boxes; returns 0 on success
int snapshot_write(const char *path, int step, double time, double dx, double dy, struct Box *const *boxes,
                   int nboxes, const int *origin_i, const int *origin_j, const int *fields, int nfields);

struct Snapshot
{
    int version, step;
    double time, dx, dy;
    int nfields, fields[SNAP_MAX_FIELDS];
    int nboxes;
    int nx[SNAP_MAX_BOXES], ny[SNAP_MAX_BOXES];
    int origin_i[SNAP_MAX_BOXES], origin_j[SNAP_MAX_BOXES];
    size_t offset[SNAP_MAX_BOXES]; // start of each box's data in data[]
    float *data;
    size_t capacity; // floats allocated in data (reused across reads)
};

// Read a snapshot into s (zero-initialize s before the first call); returns 0 on success
int snapshot_read(const char *path, struct Snapshot *s);
void snapshot_free(struct Snapshot *s);

// Index of a field in the snapshot, or -1 if it was not written
int snapshot_field_slot(const struct Snapshot *s, int field);

// Interior values of one field of one box (nx * ny, i fastest)
static inline const float *snapshot_box_field(const struct Snapshot *s, int box, int slot)
{
    return s->data + s->offset[box] + (size_t)slot * s->nx[box] * s->ny[box];
}

// Legacy text output
void save_box_to_csv(struct Box *box, const char *filename);
void save_box_connections(struct Box *boxes[], int count, const char *filename);

#endif
