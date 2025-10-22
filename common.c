#include "common.h"

uint64_t tic_time_ns = 0;

void tic(void)
{
    tic_time_ns = ns_now_monotonic();
}
time_t toc(void)
{
    return ns_now_monotonic() - tic_time_ns;
}

unsigned int verbosity_level = 0;

void printlog(unsigned int v, const char *format, ...)
{
    va_list argptr;
    if (v>verbosity_level) {
        return;
    }
    fprintf(stderr, "[%.3f] ", toc_sec());
    va_start(argptr, format);
    vfprintf(stderr, format, argptr);
    fprintf(stderr, "\n");
    va_end(argptr);
}

void increase_verbosity(void)
{
    ++verbosity_level;
}