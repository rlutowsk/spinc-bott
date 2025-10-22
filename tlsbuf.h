#ifndef TLSBUF_H
#define TLSBUF_H

#include <stddef.h>
#include <stdio.h>

/**
 * Thread-local output buffer to reduce contention on stdout.
 */
typedef struct {
    char **lines;
    size_t count;
    size_t capacity;
    FILE *out;
} OutputBuffer;

void buffer_init(OutputBuffer *buffer, size_t initial_capacity, FILE *out);

void buffer_destroy(OutputBuffer *buffer);

void buffer_add(OutputBuffer *buffer, char* line);

void buffer_flush(OutputBuffer *buffer);

#endif /* TLSBUF_H */