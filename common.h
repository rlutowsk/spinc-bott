#ifndef COMMON_H
#define COMMON_H

#include <time.h>
#include <stdio.h>

void   tic(void);
time_t toc(void);

static inline float toc_sec(void)
{
    return toc()/1000000000.0;
}

void help(const char *name);

extern unsigned int verbosity_level;

//#define log(v, ...) printf("[%.6f] ", toc_sec()); if (v>=verbosity_level) {printf(__VA_ARGS__);}

#include <stdarg.h>
void printlog(unsigned int v, const char *format, ...);

#endif