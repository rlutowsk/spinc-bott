#ifndef COMMON_H
#define COMMON_H

#define _GNU_SOURCE
#include <time.h>
#include <stdio.h>
#include <stdint.h>

#include "config.h"

#if defined(__GNUC__) || defined(__clang__)
#   define INLINE __attribute__((always_inline)) inline
#   define LIKELY(x)   __builtin_expect(!!(x), 1)
#   define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#   define INLINE inline
#   define LIKELY(x)   (x)
#   define UNLIKELY(x) (x)
#endif

#if defined(PROFILE) || defined(DEBUG)
#   ifdef INLINE
#       undef INLINE
#   endif
#   define INLINE __attribute__((noinline))
#endif

#ifndef HAVE_UINT128
#   if defined(__SIZEOF_INT128__)
#       define HAVE_UINT128 1
#   else
#       define HAVE_UINT128 0
#   endif
#endif

/*
 * set the default line length
 */
#define MAXLINE 64

/*
 * use nauty.h to decide whether we use tls or not
 */
#include <nauty.h>
#if defined(HAVE_TLS) && (HAVE_TLS == 1)
#   define NAUTY_HAS_TLS 1
#else
#   define NAUTY_HAS_TLS 0
#endif

#if NAUTY_HAS_TLS

#   include <stdatomic.h>
#   define ATOMIC_ATTR _Atomic

#   define ATOMIC_GET(x) atomic_load_explicit(&(x), memory_order_relaxed)
#   define ATOMIC_INC(x) atomic_fetch_add_explicit(&(x), 1, memory_order_relaxed)
#   define ATOMIC_DEC(x) atomic_fetch_sub_explicit(&(x), 1, memory_order_relaxed)

#else

#   define ATOMIC_ATTR

#   define ATOMIC_GET(x) (x)
#   define ATOMIC_INC(x) ((x)++)
#   define ATOMIC_DEC(x) ((x)--)

#endif

static INLINE uint64_t ns_now_monotonic(void) {
    struct timespec ts;
#if defined(HAVE_CLOCK_MONOTONIC_RAW)
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

static INLINE float toc_sec(void)
{
    return toc()/1000000000.0;
}

extern unsigned int verbosity_level;

#include <stdarg.h>
void printlog(unsigned int v, const char *format, ...);

void increase_verbosity(void);

static INLINE void remove_newline(char *s)
{
    s[strcspn(s, "\r\n")] = 0;
}

#endif