#include "common.h"

struct timespec tv;
void tic(void)
{
    clock_gettime(CLOCK_REALTIME, &tv);
}
time_t toc(void)
{
    static struct timespec st;
    clock_gettime(CLOCK_REALTIME, &st);
    return (tv.tv_sec==st.tv_sec) ? st.tv_nsec - tv.tv_nsec : 1000000000*(st.tv_sec-tv.tv_sec)+st.tv_nsec - tv.tv_nsec;
}

unsigned int verbosity_level = 0;

void printlog(unsigned int v, const char *format, ...)
{
    va_list argptr;
    if (v>verbosity_level) {
        return;
    }
    printf("[%.6f] ", toc_sec());
    va_start(argptr, format);
    vfprintf(stdout, format, argptr);
    va_end(argptr);
}