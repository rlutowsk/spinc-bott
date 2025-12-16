#pragma once

#define GINT_TO_POINTER(x) ((void*)(int)(x))

#include <stdbool.h>
#include <math.h>
#include <xxhash.h>

#include "common.h"

/* ---- Internal flat open-addressing set for 16/32-byte keys ---- */
typedef struct {
    key128_t *keys;   /* inline keys */
    uint8_t  *ctrl;   /* 0=empty, 1=full, 2=deleted */
    size_t    cap;    /* power of two */
    size_t    size;   /* # FULL slots */
    size_t    dels;   /* # tombstones */
} FlatSet;

void flat_init(FlatSet *s, size_t cap_hint);
void flat_free(FlatSet *s);
bool flat_lookup(const FlatSet *s, const key128_t *k);
bool flat_remove(FlatSet *s, const key128_t *k);
bool flat_insert(FlatSet *s, const key128_t *k);

/* Opaque bucket type with two backends:
 *  - FLAT128 (open addressing; compact) when used with key128_hash/key128_equal
 */
typedef struct _GHashBucket GHashBucket;

/* Create a sharded hash "bucket".
 * If (hash_func == key128_hash && key_equal_func == key128_equal),
 * the bucket uses FLAT128 (compact) per shard; otherwise GLib tables.
 *
 * key_destroy_func:
 *  - FLAT128 mode: used only by g_bucket_insert (pointer variant) to free caller key
 *    immediately after copying bytes internally (so callers may pass malloc'd keys).
 */

typedef void (*destroy_func_t)(void *pointer);

GHashBucket* g_bucket_new_128(
    destroy_func_t key_destroy_func,
    size_t         shards
);

/* Destroy the whole bucket and free internal storage. */
void    g_bucket_destroy (GHashBucket *b);

/* Insert (pointer API). Returns TRUE iff inserted (key was not present).
 * In FLAT128 mode, the key bytes are copied; if key_destroy_func != NULL the
 * caller pointer is freed immediately.
 */
bool g_bucket_insert (GHashBucket *b, const key128_t* key, void *value);

/* By-value insert for 16-byte keys (recommended for performance & no allocs).
 * Works in both modes:
 *  - FLAT128: copies bytes inline (no allocations).
 */
bool g_bucket_insert_copy128 (GHashBucket *b, const key128_t *key);

/* Remove key (match by value). Returns TRUE iff the key was present and removed. */
bool g_bucket_remove (GHashBucket *b, const key128_t* key);

/* Lookup key. Returns non-NULL iff present (compatible with previous usage). */
void *g_bucket_lookup (GHashBucket *b, const key128_t* key);

/* Number of elements (sum over shards). */
size_t    g_bucket_size   (GHashBucket *b);

/* Iterate all keys (16-byte) — only meaningful for FLAT128 mode.
 * 'func' is called with the address of each inline key.
 * Safe under concurrent use because we lock shards internally for iteration.
 */
typedef void (*GKey128ForeachFunc)(const key128_t *key, void *user_data);
void     g_bucket_foreach128 (GHashBucket *b, GKey128ForeachFunc func, void *user_data);

/* Hash & equal for 16-byte keys. */
static INLINE unsigned int key128_hash(const unsigned char *p) {
    return (unsigned int)XXH3_64bits(p, 16);
}

static INLINE bool key128_equal(const key128_t* a, const key128_t* b) {
#if HAVE_UINT128
    return a->u == b->u;
#else
    return memcmp(a, b, 16) == 0;
#endif
}

static INLINE bool key128_lt(const key128_t* a, const key128_t* b) {
#if HAVE_UINT128
    return a->u < b->u;
#else
    return memcmp(a, b, 16) < 0;
#endif
}

/* Mode probe (optional; handy for assertions in debugging). */
bool g_bucket_is_flat128 (GHashBucket *b);

void g_bucket_reserve(GHashBucket *b, uint64_t total_expected);

static INLINE size_t shards_for_cap_limit(uint64_t total_expected,
                     double   skew,        // np. 1.15
                     double   max_true,    // np. 0.80
                     unsigned cap_limit_pow2)
{
    // unsigned cap_limit_pow2 = pick_k_balanced(total_expected, skew, max_true, 17, 21, 128, 1, 21);
    // printlog(1, "k = %u", cap_limit_pow2);
    const double cap_limit = (double)(1u << cap_limit_pow2); // 2^k
    // Minimalna liczba shardów: ceil( total * skew / (cap_limit * max_true) )
    double need = (skew * (double)total_expected) / (cap_limit * max_true);
    size_t shards = (size_t)ceil(need);
    return shards>127 ? shards : 823;
}
