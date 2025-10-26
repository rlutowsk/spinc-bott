#ifndef TLSBUF_H
#define TLSBUF_H

#include "common.h"

/**
 * Thread-local output buffer to reduce contention on stdout.
 */
typedef struct _OutputBuffer {
    char **lines;
    size_t count;
    size_t capacity;
    FILE *out;
} OutputBuffer;

void buffer_init(OutputBuffer *buffer, size_t initial_capacity, FILE *out);

void buffer_destroy(OutputBuffer *buffer);

void buffer_add(OutputBuffer *buffer, char* line);

void buffer_flush(OutputBuffer *buffer);

void monitor_progress(unsigned long *progress_ptr, unsigned long total, FILE *stream);

#endif /* TLSBUF_H */