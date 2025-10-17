// d6pack11.h (with n stored in key)
// Ultra-fast (n <= 11) digraph6 <-> 16-byte key packing, with n embedded in key128_t
// License: MIT
//
// Layout (canonical):
//   key[127:124] : 4-bit n (1..11)
//   key[L-1:0]   : L = n*n bits of adjacency bitstream (row-major, like digraph6 R(x))
//   key[123:L]   : zero padding
//
// Reference: Brendan McKay, formats.txt: "&" + N(n) + R(x), R(x) packs row-major n*n bits in 6-bit chunks.
// https://users.cecs.anu.edu.au/~bdm/data/formats.txt

#ifndef D6PACK11_H
#define D6PACK11_H

#include <stdint.h>
#include <string.h>
#include <assert.h>

#if defined(__GNUC__) || defined(__clang__)
#  define D6PACK_INLINE __attribute__((always_inline)) inline
#  define D6PACK_LIKELY(x)   __builtin_expect(!!(x), 1)
#  define D6PACK_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#  define D6PACK_INLINE inline
#  define D6PACK_LIKELY(x)   (x)
#  define D6PACK_UNLIKELY(x) (x)
#endif

#if defined(PROFILE) && defined(D6PACK_INLINE)
#  undef D6PACK_INLINE
#  define D6PACK_INLINE __attribute__((noinline))
#endif

#ifndef D6PACK_HAVE_UINT128
#  if defined(__SIZEOF_INT128__)
#    define D6PACK_HAVE_UINT128 1
#  else
#    define D6PACK_HAVE_UINT128 0
#  endif
#endif

#ifndef D6PACK_KEY128_T_DEFINED
#define D6PACK_KEY128_T_DEFINED 1
typedef union { 
    unsigned char b[16];
#if defined(__SIZEOF_INT128__)
    __uint128_t   u;
#endif
} key128_t;
#endif

// n nibble location (byte-agnostic, portable)
#define D6PACK_N_BYTE   15u
#define D6PACK_N_SHIFT  4u
#define D6PACK_N_MASK   0x0Fu

static D6PACK_INLINE unsigned d6pack_get_n(const key128_t *k)
{ return (unsigned)((k->b[D6PACK_N_BYTE] >> D6PACK_N_SHIFT) & D6PACK_N_MASK); }

static D6PACK_INLINE void d6pack_set_n(key128_t *k, unsigned n)
{
    k->b[D6PACK_N_BYTE] =
        (unsigned char)((k->b[D6PACK_N_BYTE] & 0x0F) | (unsigned char)((n & D6PACK_N_MASK) << D6PACK_N_SHIFT));
}

// Zero bits above n^2 (keeps n in top nibble); produces canonical key
static D6PACK_INLINE void d6pack_zero_above_nsquare(key128_t *k, unsigned n)
{
    unsigned L = n * n;
    if (L >= 128) return; // for n<=11, never happens
    uint64_t lo, hi; memcpy(&lo, k->b + 0, 8); memcpy(&hi, k->b + 8, 8);
    if (L <= 64) {
        uint64_t mask = (L==64)?~0ull:((1ull << L) - 1ull);
        lo &= mask;
        hi = 0; // we'll restore n nibble below
    } else {
        unsigned r = L - 64;
        uint64_t mask = (r==64)?~0ull:((1ull << r) - 1ull);
        hi &= mask;
    }
    memcpy(k->b + 0, &lo, 8); memcpy(k->b + 8, &hi, 8);
    // restore n nibble (in case we cleared b[15])
    // unsigned n0 = d6pack_get_n(k);
    // d6pack_set_n(k, n0 ? n0 : n);
}

static D6PACK_INLINE unsigned d6pack_expected_len(unsigned n)
{
    unsigned L = n * n; unsigned groups = (L + 5) / 6;
    return 2u + groups; // '&' + N(n) + groups
}

// ---- Decode: digraph6 -> key (stores n in top nibble); returns 1/0. ----
static D6PACK_INLINE int d6pack_decode(const char *d6, key128_t *out, unsigned *out_n)
{
    if (D6PACK_UNLIKELY(!d6 || !out)) return 0;
    const unsigned char *p = (const unsigned char*)d6;
    if (D6PACK_UNLIKELY(p[0] != '&')) return 0;
    unsigned c1 = p[1]; if (D6PACK_UNLIKELY(c1 < 63)) return 0;
    unsigned n = c1 - 63u; if (D6PACK_UNLIKELY(n == 0 || n > 11)) return 0;

    unsigned L = n * n, groups = (L + 5) / 6, last = L % 6;
    unsigned need = 2 + groups;
    // hard length-check: input must be exactly 'need' chars
    for (unsigned i=0;i<need;i++) if (!p[i]) return 0;
    if (p[need] != 0) return 0; // allow only properly NUL-terminated fixed-length strings

#if D6PACK_HAVE_UINT128
    __uint128_t acc = 0;
    for (unsigned i=0;i<groups;++i) {
        unsigned v = (unsigned)p[2+i] - 63u; // 0..63
        unsigned t = (i==groups-1 && last) ? last : 6u;
        acc = (acc << t) | (unsigned)(v >> (6u - t));
    }
    // uint64_t lo = (uint64_t)acc, hi = (uint64_t)(acc >> 64);
    // memcpy(out->b+0,&lo,8); memcpy(out->b+8,&hi,8);
    out->u = acc;
#else
    uint64_t lo=0, hi=0; unsigned pos=0;
    for (unsigned i=0;i<groups && pos<L;++i) {
        unsigned v = (unsigned)p[2+i] - 63u; unsigned t = (i==groups-1 && last) ? last : 6u;
        unsigned bits = v >> (6u - t); unsigned left = t;
        while (left) {
            unsigned lane_rem = (pos<64)? (64u-pos) : (128u-pos);
            unsigned take     = (left < lane_rem)? left : lane_rem;
            uint64_t chunk = (uint64_t)(bits >> (left - take)) & ((take==64)?~0ull:((1ull<<take)-1ull));
            if (pos<64) lo |= chunk << pos; else hi |= chunk << (pos-64);
            pos += take; left -= take;
        }
    }
    memcpy(out->b+0,&lo,8); memcpy(out->b+8,&hi,8);
#endif
    d6pack_zero_above_nsquare(out, n); // canon: zero padding
    d6pack_set_n(out, n);              // put n nibble in b[15]
    if (out_n) *out_n = n;
    return 1;
}

// ---- Encode from key (uses n stored in nibble). Returns chars without '\0'. ----
static D6PACK_INLINE unsigned d6pack_encode_from_key(const key128_t *k, char *out_d6)
{
    unsigned n = d6pack_get_n(k);
    if (D6PACK_UNLIKELY(n < 1 || n > 11)) return 0;
    unsigned L = n * n, last = L % 6;

    char *q = out_d6; *q++='&'; *q++=(char)(63+n);
#if D6PACK_HAVE_UINT128
    uint64_t lo, hi; memcpy(&lo, k->b+0, 8); memcpy(&hi, k->b+8, 8);
    __uint128_t acc = ((__uint128_t)hi << 64) | lo;
    int bitpos = (int)L;
    for (unsigned i = 0; i < (L/6); ++i) { bitpos -= 6; unsigned v = (unsigned)((acc >> bitpos) & 0x3Fu); *q++ = (char)(63 + v); }
    if (last) { bitpos -= (int)last; unsigned v = (unsigned)((acc >> bitpos) & ((1u << last) - 1u)); v <<= (6u - last); *q++ = (char)(63 + v); }
#else
    uint64_t lo, hi; memcpy(&lo, k->b+0, 8); memcpy(&hi, k->b+8, 8);
    for (unsigned g = 0; g < (L/6); ++g) {
        unsigned base = L - 6*(g+1); unsigned v = 0;
        for (unsigned j=0;j<6;++j) {
            unsigned idx = base + (5-j);
            unsigned b = (idx<64u)? ((lo>>idx)&1u) : ((hi>>(idx-64u))&1u);
            v = (v<<1) | b;
        }
        *q++ = (char)(63 + (v & 0x3F));
    }
    if (last) {
        unsigned v = 0; for (unsigned j=0;j<last;++j) { unsigned b = (j<64u)? ((lo>>j)&1u) : ((hi>>(j-64u))&1u); v |= (b << (last-1-j)); }
        v <<= (6 - last); *q++ = (char)(63 + (v & 0x3F));
    }
#endif
    *q = '\0';
    return (unsigned)(q - out_d6);
}

// Convenience: encode with explicit n (nie modyfikuje oryginalnego klucza)
static D6PACK_INLINE unsigned d6pack_encode(const key128_t *k, unsigned n, char *out_d6)
{
    key128_t tmp = *k; d6pack_set_n(&tmp, n); d6pack_zero_above_nsquare(&tmp, n);
    return d6pack_encode_from_key(&tmp, out_d6);
}

#endif // D6PACK11_H