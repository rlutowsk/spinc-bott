#include <omp.h>

#include "tlsbuf.h"

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

/* --- Progress monitor: prints to the specified stream (stderr with -o) --- */
void monitor_progress(unsigned long *progress_ptr, unsigned long total, FILE *stream)
{
    const double period = 0.25;
    double next = omp_get_wtime();
    double start= omp_get_wtime();
    double old  = start;
    struct timespec ts = { .tv_sec = 0, .tv_nsec = 500000000 }; /* 500 ms */
    unsigned long cur;
    double now, pct, speed;
    char etabuf[256] = {0};

    for (;;) {
        #pragma omp atomic read
            cur = *progress_ptr;

        if (cur >= total) break;

        now = omp_get_wtime();
        if (now >= next) {
            pct   = (total == 0) ? 100.0 : (100.0 * (double)cur / (double)total);
            if (now >= old + 2.5) {
                old = now;
                speed = pct / (now-start);
                unsigned long s = (unsigned long)((100.0-pct)/speed);
                unsigned long h = s / 3600;
                unsigned long m = ( s % 3600) / 60;
                snprintf(etabuf, 255, " (ETA %lu:%02lu:%02lu)", h, m, s % 60);
            }
            fprintf(stream, "\33[2K\r%8.4f%% %s", pct, etabuf);
            fflush(stream);
            next = now + period;
        }
        nanosleep(&ts, NULL);
    }
}