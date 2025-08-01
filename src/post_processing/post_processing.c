#include <stdio.h>
#include <stdlib.h>

#include "post_processing.h"

void save_box_to_csv(struct Box *box, const char *filename)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("Error opening file");
        return;
    }

    int nx_total = box->nx + 2;

    // CSV header
    fprintf(fp, "i,j,w_1,w_2,w_3,w_4,u_1,u_2,u_3,u_4\n");

    // Loop over interior cells only (skip ghost cells)
    for (int j = 1; j <= box->ny; ++j)
    {
        for (int i = 1; i <= box->nx; ++i)
        {
            int idx = j * nx_total + i;

            fprintf(fp, "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n",
                    i - 1, j - 1, // shift index to start from 0
                    box->w_1[idx], box->w_2[idx], box->w_3[idx], box->w_4[idx],
                    box->u_1[idx], box->u_2[idx], box->u_3[idx], box->u_4[idx]);
        }
    }

    fclose(fp);
}

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

            const char *type_str = (bc->type == BOUNDARY_REFLECTIVE) ? "REFLECTIVE" : (bc->type == BOUNDARY_TRANSMISSIVE) ? "TRANSMISSIVE"
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
