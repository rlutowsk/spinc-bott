#ifndef BUCKET_H
#define BUCKET_H

#include <glib.h>
#include <stdbool.h>
#include <math.h>
#include <xxhash.h>

#if HASH_MIXED_MODE == 1
#   define __USE_GLIB
#else
#   undef __USE_GLIB
#endif

#include "d6pack11.h"
/* ---- 16-byte key type for compact exact IDs (canonical d6 zero-padded) ---- */
#ifndef K_LEN
#   define K_LEN 16
#endif
#ifndef D6PACK_KEY128_T_DEFINED
typedef union { 
    unsigned char b[K_LEN];
#if defined(__SIZEOF_INT128__)
    __uint128_t   u;
#endif
} key128_t;
#endif

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
 *  - GLIB (GHashTable shards) otherwise
 */
typedef struct _GHashBucket GHashBucket;

/* Create a sharded hash "bucket".
 * If (hash_func == key128_hash && key_equal_func == key128_equal),
 * the bucket uses FLAT128 (compact) per shard; otherwise GLib tables.
 *
 * key_destroy_func:
 *  - GLIB mode: stored and used by GHashTable to free keys on remove/destroy.
 *  - FLAT128 mode: used only by g_bucket_insert (pointer variant) to free caller key
 *    immediately after copying bytes internally (so callers may pass malloc'd keys).
 * value_destroy_func: only used in GLIB mode (we do not store values in FLAT128).
 */
// GHashBucket* g_bucket_new_full(
//     GHashFunc      hash_func,
//     GEqualFunc     key_equal_func,
//     GDestroyNotify key_destroy_func,
//     GDestroyNotify value_destroy_func,
//     gsize          shards
// );

GHashBucket* g_bucket_new_128(
    GDestroyNotify key_destroy_func,
    GDestroyNotify value_destroy_func,
    gsize          shards
);

/* Destroy the whole bucket and free internal storage. */
void    g_bucket_destroy (GHashBucket *b);

/* Insert (pointer API). Returns TRUE iff inserted (key was not present).
 * In FLAT128 mode, the key bytes are copied; if key_destroy_func != NULL the
 * caller pointer is freed immediately. In GLIB mode the key pointer becomes
 * owned by the table (freed via key_destroy_func on removal/destroy).
 */
gboolean g_bucket_insert (GHashBucket *b, const key128_t* key, gpointer value);

/* By-value insert for 16-byte keys (recommended for performance & no allocs).
 * Works in both modes:
 *  - FLAT128: copies bytes inline (no allocations).
 *  - GLIB: duplicates a 16-byte heap key and inserts with (value == TRUE).
 */
gboolean g_bucket_insert_copy128 (GHashBucket *b, const key128_t *key);

/* Remove key (match by value). Returns TRUE iff the key was present and removed. */
gboolean g_bucket_remove (GHashBucket *b, const key128_t* key);

/* Lookup key. Returns non-NULL iff present (compatible with previous usage). */
gpointer g_bucket_lookup (GHashBucket *b, const key128_t* key);

/* Number of elements (sum over shards). */
gsize    g_bucket_size   (GHashBucket *b);

/* Iterate all keys (16-byte) — only meaningful for FLAT128 mode.
 * 'func' is called with the address of each inline key.
 * Safe under concurrent use because we lock shards internally for iteration.
 */
typedef void (*GKey128ForeachFunc)(const key128_t *key, gpointer user_data);
void     g_bucket_foreach128 (GHashBucket *b, GKey128ForeachFunc func, gpointer user_data);

/* Hash & equal for 16-byte keys. */
static inline guint key128_hash(const unsigned char *p) {
    return (guint)XXH3_64bits(p, 16);
}

static inline gboolean key128_equal(const key128_t* a, const key128_t* b) {
#if defined(__SIZEOF_INT128__)
    return a->u == b->u;
#else
    return memcmp(a, b, 16) == 0;
#endif
}

static inline gboolean key128_lt(const key128_t* a, const key128_t* b) {
#if defined(__SIZEOF_INT128__)
    return a->u < b->u;
#else
    return memcmp(a, b, 16) < 0;
#endif
}

/* Mode probe (optional; handy for assertions in debugging). */
gboolean g_bucket_is_flat128 (GHashBucket *b);

void g_bucket_reserve(GHashBucket *b, uint64_t total_expected);

static inline size_t shards_for_cap_limit(uint64_t total_expected,
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

#endif /* BUCKET_H */
