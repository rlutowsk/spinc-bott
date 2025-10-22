#ifndef COMMON_H
#define COMMON_H

#define _GNU_SOURCE
#include <time.h>
#include <stdio.h>
#include <stdint.h>

static inline uint64_t ns_now_monotonic(void) {
    struct timespec ts;
#if defined(CLOCK_MONOTONIC_RAW)
    const clockid_t cid = CLOCK_MONOTONIC_RAW;
#else
    const clockid_t cid = CLOCK_MONOTONIC;
#endif    
    if (clock_gettime(cid, &ts) != 0) {
            return 0;
    }
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}

void   tic(void);
time_t toc(void);

static inline float toc_sec(void)
{
    return toc()/1000000000.0;
}

extern unsigned int verbosity_level;

#include <stdarg.h>
void printlog(unsigned int v, const char *format, ...);

void increase_verbosity(void);

#endif