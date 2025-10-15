#ifndef ADJPACK11_H
#define ADJPACK11_H

#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "d6pack11.h"
#include "bott.h"

#if defined(__GNUC__) || defined(__clang__)
# define ADJPACK_INLINE __attribute__((always_inline)) inline
# define ADJPACK_LIKELY(x) __builtin_expect(!!(x), 1)
# define ADJPACK_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
# define ADJPACK_INLINE inline
# define ADJPACK_LIKELY(x) (x)
# define ADJPACK_UNLIKELY(x) (x)
#endif

#define ADJPACK_N_BYTE 15u
#define ADJPACK_N_SHIFT 4u
#define ADJPACK_N_MASK 0x0Fu

static ADJPACK_INLINE void adjpack_set_n(key128_t *k, unsigned n) {
    k->b[ADJPACK_N_BYTE] =
        (unsigned char)((k->b[ADJPACK_N_BYTE] & 0x0F) |
        ((n & ADJPACK_N_MASK) << ADJPACK_N_SHIFT));
}

static ADJPACK_INLINE unsigned adjpack_get_n(const key128_t *k) {
    return (unsigned)((k->b[ADJPACK_N_BYTE] >> ADJPACK_N_SHIFT) & ADJPACK_N_MASK);
}

static ADJPACK_INLINE void adjpack_zero_above_nsquare(key128_t *k, unsigned n) {
    unsigned L = n * n;
    if (L >= 128) return;
#if D6PACK_HAVE_UINT128
	k->u &= (((__uint128_t)1 << L) - 1);
#else
    uint64_t lo, hi;
    memcpy(&lo, k->b + 0, 8);
    memcpy(&hi, k->b + 8, 8);
    if (L <= 64) {
        uint64_t mask = (L == 64) ? ~0ull : ((1ull << L) - 1ull);
        lo &= mask;
        hi = 0;
    } else {
        unsigned r = L - 64;
        uint64_t mask = (r == 64) ? ~0ull : ((1ull << r) - 1ull);
        hi &= mask;
    }
    memcpy(k->b + 0, &lo, 8);
    memcpy(k->b + 8, &hi, 8);
#endif
}

static ADJPACK_INLINE uint64_t bitreverse64(uint64_t x) {
    x = ((x & 0x5555555555555555ULL) <<  1) | ((x >>  1) & 0x5555555555555555ULL);
    x = ((x & 0x3333333333333333ULL) <<  2) | ((x >>  2) & 0x3333333333333333ULL);
    x = ((x & 0x0F0F0F0F0F0F0F0FULL) <<  4) | ((x >>  4) & 0x0F0F0F0F0F0F0F0FULL);
    x = ((x & 0x00FF00FF00FF00FFULL) <<  8) | ((x >>  8) & 0x00FF00FF00FF00FFULL);
    x = ((x & 0x0000FFFF0000FFFFULL) << 16) | ((x >> 16) & 0x0000FFFF0000FFFFULL);
    x = (x << 32) | (x >> 32);
    return x;
}

// Packs adjacency matrix mat[n] into key128_t (compatible with digraph6)
static ADJPACK_INLINE void adjpack_from_matrix(const vec_t *mat, unsigned n, key128_t *out) {
    assert(n >= 1 && n <= 11);
#if D6PACK_HAVE_UINT128
	out->u = 0;
    for (unsigned i = 0; i < n; ++i) {
        uint64_t row = bitreverse64(mat[i]) >> (64 - n);
        out->u <<= n;
        out->u |= row;
    }
#else
    uint64_t lo = 0, hi = 0;
    unsigned pos = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (int j = n - 1; j >= 0; --j) {
            uint64_t bit = (mat[i] >> j) & 1;
            if (pos < 64)
                lo |= bit << pos;
            else
                hi |= bit << (pos - 64);
            ++pos;
        }
    }
    memcpy(out->b + 0, &lo, 8);
    memcpy(out->b + 8, &hi, 8);
#endif
    adjpack_zero_above_nsquare(out, n);
    adjpack_set_n(out, n);
}


// Unpacks key128_t into adjacency matrix mat[n]
static ADJPACK_INLINE void adjpack_to_matrix(const key128_t *k, uint64_t *mat_out, unsigned n) {
    assert(n >= 1 && n <= 11);
    //memset(mat_out, 0, n * sizeof(uint64_t));
#if D6PACK_HAVE_UINT128
    /*
    __uint128_t acc = k->u;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            unsigned bitpos = n * n - 1 - (i * n + j);
            uint64_t bit = (acc >> bitpos) & 1;
            mat_out[i] |= bit << j;
        }
    }
    */
    __uint128_t acc = k->u;
    for (int i = n-1; i >= 0; --i) {
        mat_out[i] = bitreverse64(acc) >> (64 - n);
        acc >>= n;
    }
#else
    uint64_t lo, hi;
    memcpy(&lo, k->b + 0, 8);
    memcpy(&hi, k->b + 8, 8);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            unsigned bitpos = i * n + j;
            uint64_t bit = (bitpos < 64) ? ((lo >> bitpos) & 1) : ((hi >> (bitpos - 64)) & 1);
            mat_out[i] |= bit << j;
        }
    }
#endif
}

#endif // ADJPACK11_H
