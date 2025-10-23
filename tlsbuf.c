#include "tlsbuf.h"

#include <stdlib.h>
#include <string.h>
#include <omp.h>

void buffer_init(OutputBuffer *buffer, size_t initial_capacity, FILE *out) {
    buffer->lines = malloc(initial_capacity * sizeof(char*));
    buffer->count = 0;
    buffer->capacity = initial_capacity;
    buffer->out = out;
}
void buffer_destroy(OutputBuffer *buffer) {
    // if flush() has been called, count == 0 and this is a no-op
    // for (size_t i = 0; i < buffer->count; ++i) free(buffer->lines[i]);
    free(buffer->lines);
    buffer->lines = NULL; buffer->count = 0; buffer->capacity = 0;
}
void buffer_add(OutputBuffer *buffer, char* line) {
    if (buffer->count >= buffer->capacity) {
        #pragma omp critical (StdoutLock)
        {
            for (size_t i = 0; i < buffer->count; ++i) {
                fputs(buffer->lines[i], buffer->out);
                if (buffer->lines[i][strlen(buffer->lines[i])-1] != '\n')
                    fputc('\n', buffer->out);
                //free(buffer->lines[i]);
            }
        }
        buffer->count = 0;
    }
    buffer->lines[buffer->count++] = line; //strdup(line);
}
void buffer_flush(OutputBuffer *buffer) {
    if (buffer->count > 0) {
        #pragma omp critical (StdoutLock)
        {
            for (size_t i = 0; i < buffer->count; ++i) {
                fputs(buffer->lines[i], buffer->out);
                if (buffer->lines[i][strlen(buffer->lines[i])-1] != '\n')
                    fputc('\n', buffer->out);
                //free(buffer->lines[i]);
            }
        }
    }
    buffer->count = 0;
}