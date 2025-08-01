#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H

#include "../box/box.h"

void save_box_to_csv(struct Box *box, const char *filename);

void save_box_connections(struct Box *boxes[], int count, const char *filename);

#endif
