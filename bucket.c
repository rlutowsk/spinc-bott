#include "common.h"

#include <omp.h>

#if NAUTY_HAS_TLS

#   define MUTEX_LOCK(x)   omp_set_lock(x)
#   define MUTEX_UNLOCK(x) omp_unset_lock(x)

#else

#   define MUTEX_LOCK(x)
#   define MUTEX_UNLOCK(x)

#endif

#include "bucket.h"

#ifndef MIN_TRUE
#define MIN_TRUE 0.20
#endif

#ifndef MAX_TRUE
#define MAX_TRUE 0.80
#endif

#ifndef MAX_OCCUPIED
#define MAX_OCCUPIED 0.80
#endif

#define CTRL_EMPTY   ((uint8_t)0)
#define CTRL_FULL    ((uint8_t)1)
#define CTRL_DELETED ((uint8_t)2)

static INLINE uint64_t hash16(const key128_t *k)
{
    return XXH3_64bits(k->b, 16);
}

static size_t next_pow2(size_t x) {
    // Keep the original behavior for small values.
    if (x <= 8) return 8;
    x--;
#if SIZE_MAX > 0xFFFFFFFFu
    // 64-bit system: use __builtin_clzl (for unsigned long)
    return (size_t)1 << (64 - __builtin_clzl(x));
#else
    // 32-bit system: use __builtin_clz (for unsigned int)
    return (size_t)1 << (32 - __builtin_clz(x));
#endif
}

void flat_init(FlatSet *s, size_t cap_hint) {
    s->cap  = next_pow2(cap_hint ? cap_hint : 1024);
    s->size = 0;
    s->dels = 0;
    s->keys = (key128_t*)calloc(s->cap, sizeof(key128_t));
    s->ctrl = (uint8_t*) calloc(s->cap, sizeof(uint8_t) );
}

void flat_free(FlatSet *s) {
    if (!s){
        return;
    }
    free(s->keys);
    free(s->ctrl);
    s->keys = NULL;
    s->ctrl = NULL;
    s->cap = s->size = s->dels = 0;
}

static void flat_rehash(FlatSet *s, size_t new_cap) {
    FlatSet dst;
    dst.cap = 0;
    dst.size = dst.dels = 0;
    dst.keys = NULL;
    dst.ctrl = NULL;
    flat_init(&dst, new_cap);

    for (size_t i = 0; i < s->cap; ++i) {
        if (s->ctrl[i] == CTRL_FULL) {
            key128_t k = s->keys[i];
            uint64_t h = hash16(&k);
            size_t m = dst.cap - 1;
            size_t j = (size_t)(h & m);
            while (dst.ctrl[j] == CTRL_FULL) {
                j = (j + 1) & m;
            }
            dst.keys[j] = k;
            dst.ctrl[j] = CTRL_FULL;
            dst.size++;
        }
    }
    free(s->keys);
    free(s->ctrl);
    *s = dst;
}

bool flat_lookup(const FlatSet *s, const key128_t *k) {
    if (s->cap == 0) return false;
    uint64_t h = hash16(k);
    size_t m = s->cap - 1, i = (size_t)(h & m);
    for (;;) {
        uint8_t c = s->ctrl[i];
        if (c == CTRL_EMPTY) {
            return false;
        }
        if (c == CTRL_FULL && key128_equal(&s->keys[i], k)){
            return TRUE;
        }
        i = (i + 1) & m;
    }
}


static INLINE double flat_load(const FlatSet *s) {
    return (double)(s->size + s->dels) / (double)s->cap;
}

// --- New: separate real (true) load and occupied load (FULL + DELETED)
static INLINE double true_load(const FlatSet *s) {
    return s->cap ? (double)s->size / (double)s->cap : 0.0;
}
static INLINE double occupied_load(const FlatSet *s) {
    return s->cap ? (double)(s->size + s->dels) / (double)s->cap : 0.0;
}

// --- New: optional shrink with hysteresis.
// Shrink only when true load is clearly low, to avoid ping-pong.
static void flat_maybe_shrink(FlatSet *s)
{
    const size_t MIN_CAP  = 16;

    if (s->cap > MIN_CAP && true_load(s) < MIN_TRUE) {
        size_t target = (size_t)((double)s->size / MAX_TRUE) + 1;
        if (target < MIN_CAP) target = MIN_CAP;
        if (target < s->cap) {
            flat_rehash(s, target); // aligns to next_pow2
        }
    }
}

bool flat_remove(FlatSet *s, const key128_t *k) {
    if (s->cap == 0) return false;

    uint64_t h = hash16(k);
    size_t m = s->cap - 1;
    size_t i = (size_t)(h & m);

    for (;;) {
        uint8_t c = s->ctrl[i];
        if (c == CTRL_EMPTY){
            return false;
        }
        if (c == CTRL_FULL && key128_equal(&s->keys[i], k)) {
            s->ctrl[i] = CTRL_DELETED;
            s->size--;
            s->dels++;
            if (s->dels > s->size && s->dels > s->cap / 8) {
                flat_rehash(s, s->cap);
            } else {
                // possible shrink with hysteresis (only if true load is low)
                flat_maybe_shrink(s);
            }

            return TRUE;
        }
        i = (i + 1) & m;
    }
}

bool flat_insert(FlatSet *s, const key128_t *k) {
    if (s->cap == 0) flat_init(s, 1024);

    if (occupied_load(s) > MAX_OCCUPIED) {
        if (true_load(s) <= MAX_TRUE) {
            // lots of tombstones -> compact
            flat_rehash(s, s->cap);
        } else {
            // tight -> grow
            size_t need = (size_t)((double)(s->size + 1) / MAX_TRUE) + 1;
            flat_rehash(s, need);  // next_pow2 will get the right size
        }
    }

    uint64_t h = hash16(k);
    size_t m = s->cap - 1, i = (size_t)(h & m), first_del = (size_t)(-1);
    for (;;) {
        uint8_t c = s->ctrl[i];
        if (c == CTRL_EMPTY) {
            size_t pos = (first_del != (size_t)(-1)) ? first_del : i;
            s->keys[pos] = *k;
            s->ctrl[pos] = CTRL_FULL;
            s->size++;
            if (first_del != (size_t)(-1)) {
                s->dels--;
            }
            return TRUE;
        }
        if (c == CTRL_DELETED) {
            if (first_del == (size_t)(-1)) {
                first_del = i;
            }
        } else {
            if (key128_equal(&s->keys[i], k)) {
                return false;
            }
        }
        i = (i + 1) & m;
    }
}

// --- New: ensure capacity for expected elements at target max_true load.
// Rounds up to next_pow2 via flat_rehash()->flat_init()
static void flat_reserve(FlatSet *s, size_t n_expected, double max_true) {
    if (n_expected == 0) return;

    // how many slots are needed so that (size / cap) <= max_true
    size_t need = (size_t)((double)n_expected / max_true) + 1;

    if (s->cap < need) {
        flat_rehash(s, need);  // will round up to the nearest power of two
    } else if (s->dels > 0) {
        // It's worth compacting tombstones if we already have capacity
        flat_rehash(s, s->cap);
    }
}

/* ---- Bucket with shard locks ---- */
struct _GHashBucket {
    size_t      nshards;
    ATOMIC_ATTR size_t number_of_elements;

    /* Locks: one per shard (used in both modes). */
    omp_lock_t    *locks;

    /* release the elements */
    destroy_func_t key_destroy_func;

    FlatSet   *flat;
};

static INLINE size_t shard_index_flat(const GHashBucket *b, const key128_t *k) {
    return (size_t)(hash16(k) % (uint64_t)b->nshards);
}
static INLINE size_t shard_index_flat_h(const GHashBucket *b, uint64_t h) {
    return (size_t)(h % (uint64_t)b->nshards);
}

GHashBucket* g_bucket_new_128(
    destroy_func_t key_destroy_func,
    size_t         shards)
{
    GHashBucket *bucket = (GHashBucket*)malloc(sizeof(GHashBucket));
    bucket->nshards = shards ? shards : 1;
    bucket->number_of_elements = 0;
    bucket->key_destroy_func = key_destroy_func;
    bucket->locks = (omp_lock_t*)calloc(bucket->nshards, sizeof(omp_lock_t));
    for (size_t i = 0; i < bucket->nshards; ++i) omp_init_lock(&bucket->locks[i]);

    bucket->flat = (FlatSet*)calloc(bucket->nshards, sizeof(FlatSet));

    return bucket;
}

void g_bucket_destroy(GHashBucket *b)
{
    if (!b) return;
    for (size_t i = 0; i < b->nshards; ++i) flat_free(&b->flat[i]);
    free(b->flat);
    for (size_t i = 0; i < b->nshards; ++i) omp_destroy_lock(&b->locks[i]);
    free(b->locks);
    free(b);
}

size_t g_bucket_size(GHashBucket *b) {
    return b ? ATOMIC_GET(b->number_of_elements) : 0;
}

void *g_bucket_lookup(GHashBucket* b, const key128_t* k)
{
    if (!b) return NULL;
    size_t idx = shard_index_flat(b, k);
    MUTEX_LOCK(&b->locks[idx]);
    bool ok = flat_lookup(&b->flat[idx], k);
    MUTEX_UNLOCK(&b->locks[idx]);
    return ok ? GINT_TO_POINTER(TRUE) : NULL;
}

bool g_bucket_remove(GHashBucket* b, const key128_t* k)
{
    if (!b) return false;
    bool res = false;

    size_t idx = shard_index_flat(b, k);
    MUTEX_LOCK(&b->locks[idx]);
    res = flat_remove(&b->flat[idx], k);
    MUTEX_UNLOCK(&b->locks[idx]);
    if (res) {
        ATOMIC_DEC(b->number_of_elements);
    }

    return res;
}

bool g_bucket_insert(GHashBucket *b, const key128_t* kptr, void *value)
{
    (void)value; if (!b) return false;

    size_t idx = shard_index_flat(b, kptr);
    MUTEX_LOCK(&b->locks[idx]);

    bool ins = flat_insert(&b->flat[idx], kptr);
    MUTEX_UNLOCK(&b->locks[idx]);
    /* We copied bytes; respect caller ownership per contract. */

    if (ins) {
        ATOMIC_INC(b->number_of_elements);
    }
    return ins;
}

bool g_bucket_insert_copy128(GHashBucket *b, const key128_t *key)
{
    if (!b) return false;

    size_t idx = shard_index_flat(b, key);
    MUTEX_LOCK(&b->locks[idx]);
    if (b->flat[idx].cap == 0) flat_init(&b->flat[idx], 1024);
    bool ins = flat_insert(&b->flat[idx], key);
    MUTEX_UNLOCK(&b->locks[idx]);
    if (ins) ATOMIC_INC(b->number_of_elements);
    return ins;
}

void g_bucket_foreach128(GHashBucket *b, GKey128ForeachFunc func, void *user_data)
{
    if (!b || !func) return;

    for (size_t s = 0; s < b->nshards; ++s) {
        MUTEX_LOCK(&b->locks[s]);
        FlatSet *fs = &b->flat[s];
        for (size_t i = 0; i < fs->cap; ++i) {
            if (fs->ctrl[i] == CTRL_FULL) func(&fs->keys[i], user_data);
        }
        MUTEX_UNLOCK(&b->locks[s]);
    }
}

bool g_bucket_is_flat128(GHashBucket* b) { return b!=NULL; }

// --- New: reserve total capacity across shards.
// NOTE: call this before concurrent use. No locking inside.
void g_bucket_reserve(GHashBucket *b, uint64_t total_expected)
{
    if (!b || total_expected == 0) return;

    const double SKEW     = 1.15;

    // Simple per-shard expected count with skew factor
    double base = (double) total_expected / (double) b->nshards;
    size_t per_shard_expected = (size_t)(base * SKEW);

    for (size_t i = 0; i < b->nshards; ++i) {
        flat_reserve(&b->flat[i], per_shard_expected, MAX_TRUE);
    }
}
