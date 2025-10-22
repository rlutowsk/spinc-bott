#include "bucket.h"
#include <string.h>
#include <xxhash.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef MIN_TRUE
#define MIN_TRUE 0.20
#endif

#ifndef MAX_TRUE
#define MAX_TRUE 0.80
#endif

#ifndef MAX_OCCUPIED
#define MAX_OCCUPIED 0.80
#endif

typedef enum { BUCKET_GLIB = 0, BUCKET_FLAT128 = 1 } BucketKind;

#define CTRL_EMPTY   ((uint8_t)0)
#define CTRL_FULL    ((uint8_t)1)
#define CTRL_DELETED ((uint8_t)2)

static inline uint64_t hash16(const key128_t *k)
{
    return XXH3_64bits(k->b, 16);
}

static size_t next_pow2(size_t x) {
    // Zachowaj oryginalne zachowanie dla małych wartości.
    if (x <= 8) return 8;
    x--;
#if SIZE_MAX > 0xFFFFFFFFu
    // System 64-bitowy: używamy __builtin_clzl (dla unsigned long)
    return (size_t)1 << (64 - __builtin_clzl(x));
#else
    // System 32-bitowy: używamy __builtin_clz (dla unsigned int)
    return (size_t)1 << (32 - __builtin_clz(x));
#endif
}

void flat_init(FlatSet *s, size_t cap_hint) {
    s->cap  = next_pow2(cap_hint ? cap_hint : 1024);
    s->size = 0;
    s->dels = 0;
    s->keys = (key128_t*)g_malloc_n(s->cap, sizeof(key128_t));
    s->ctrl = (uint8_t*)g_malloc0_n(s->cap, sizeof(uint8_t));
}

void flat_free(FlatSet *s) {
    if (!s){
        return;
    }
    g_free(s->keys);
    g_free(s->ctrl);
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
    g_free(s->keys); 
    g_free(s->ctrl);
    *s = dst;
}

bool flat_lookup(const FlatSet *s, const key128_t *k) {
    if (s->cap == 0) return FALSE;
    uint64_t h = hash16(k);
    size_t m = s->cap - 1, i = (size_t)(h & m);
    for (;;) {
        uint8_t c = s->ctrl[i];
        if (c == CTRL_EMPTY) {
            return FALSE;
        }
        // if (c == CTRL_FULL && memcmp(&s->keys[i], k, K_LEN) == 0){
        if (c == CTRL_FULL && key128_equal(&s->keys[i], k)){
            return TRUE;
        }
        i = (i + 1) & m;
    }
}


static inline double flat_load(const FlatSet *s) {
    return (double)(s->size + s->dels) / (double)s->cap;
}

// --- New: separate real (true) load and occupied load (FULL + DELETED)
static inline double true_load(const FlatSet *s) {
    return s->cap ? (double)s->size / (double)s->cap : 0.0;
}
static inline double occupied_load(const FlatSet *s) {
    return s->cap ? (double)(s->size + s->dels) / (double)s->cap : 0.0;
}

// --- New: optional shrink with hysteresis.
// Shrink only when true load is clearly low, to avoid ping-pong.
static void flat_maybe_shrink(FlatSet *s)
{
    // const double MIN_TRUE = 0.20;  // poniżej tego próbujemy zmniejszać
    // const double MAX_TRUE = 0.80;  // po zmniejszeniu chcemy <= 80%
    const size_t MIN_CAP  = 16;

    if (s->cap > MIN_CAP && true_load(s) < MIN_TRUE) {
        size_t target = (size_t)((double)s->size / MAX_TRUE) + 1;
        if (target < MIN_CAP) target = MIN_CAP;
        if (target < s->cap) {
            flat_rehash(s, target); // wyrówna do najbliższej potęgi 2
        }
    }
}

bool flat_remove(FlatSet *s, const key128_t *k) {
    if (s->cap == 0) return FALSE;

    uint64_t h = hash16(k);
    size_t m = s->cap - 1;
    size_t i = (size_t)(h & m);

    for (;;) {
        uint8_t c = s->ctrl[i];
        if (c == CTRL_EMPTY){
            return FALSE;
        }
        // if (c == CTRL_FULL && memcmp(&s->keys[i], k, K_LEN) == 0) {
        if (c == CTRL_FULL && key128_equal(&s->keys[i], k)) {
            s->ctrl[i] = CTRL_DELETED;
            s->size--;
            s->dels++;
            if (s->dels > s->size && s->dels > s->cap / 8) {
                flat_rehash(s, s->cap);
            } else {
                // ewentualny shrink z histerezą (lekki, rzadko zadziała)
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
            // dużo grobów -> kompaktuj zamiast rosnąć
            flat_rehash(s, s->cap);
        } else {
            // faktycznie za ciasno -> zwiększ do rozmiaru dającego MAX_TRUE
            size_t need = (size_t)((double)(s->size + 1) / MAX_TRUE) + 1;
            flat_rehash(s, need);  // next_pow2 wykona pośrednie zaokrąglenie
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
            // if (memcmp(&s->keys[i], k, K_LEN) == 0) {
            if (key128_equal(&s->keys[i], k)) {
                return FALSE;
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

/* ---- Hybrid bucket with shard locks ---- */
struct _GHashBucket {
    BucketKind  kind;
    gsize       nshards;
    ATOMIC_ATTR gsize number_of_elements;

    /* Locks: one per shard (used in both modes). */
    GMutex    *locks;

    /* GLIB mode */
    GHashTable **tables;
    GDestroyNotify key_destroy_func;
    GDestroyNotify value_destroy_func;
    GHashFunc      shard_hash_func;

    /* FLAT128 mode */
    FlatSet   *flat;
};

static inline size_t shard_index_glib(const GHashBucket *b, gconstpointer key) {
    GHashFunc hf = b->shard_hash_func ? b->shard_hash_func : g_direct_hash;
    return hf(key) % b->nshards;
}
static inline size_t shard_index_flat(const GHashBucket *b, const key128_t *k) {
    return (size_t)(hash16(k) % (uint64_t)b->nshards);
}
static inline size_t shard_index_flat_h(const GHashBucket *b, uint64_t h) {
    return (size_t)(h % (uint64_t)b->nshards);
}

// GHashBucket* g_bucket_new_full(
//     GHashFunc      hash_func,
//     GEqualFunc     key_equal_func,
//     GDestroyNotify key_destroy_func,
//     GDestroyNotify value_destroy_func,
//     gsize          shards)
// {
//     GHashBucket *bucket = (GHashBucket*)g_malloc0(sizeof(GHashBucket));
//     bucket->nshards = shards ? shards : 1;
//     bucket->number_of_elements = 0;
//     bucket->key_destroy_func = key_destroy_func;
//     bucket->value_destroy_func = value_destroy_func;
//     bucket->locks = (GMutex*)g_malloc0_n(bucket->nshards, sizeof(GMutex));
//     for (gsize i = 0; i < bucket->nshards; ++i) g_mutex_init(&bucket->locks[i]);

//     gboolean is_key128 = (hash_func == key128_hash) && (key_equal_func == key128_equal);

//     if (is_key128) {
//         bucket->kind = BUCKET_FLAT128;
//         bucket->flat = (FlatSet*)g_malloc0_n(bucket->nshards, sizeof(FlatSet));
//         bucket->tables = NULL;
//         bucket->shard_hash_func = NULL;
//     } else {
// #ifdef __USE_GLIB
//         bucket->kind = BUCKET_GLIB;
//         bucket->tables = (GHashTable**)g_malloc_n(bucket->nshards, sizeof(GHashTable*));
//         for (gsize i = 0; i < bucket->nshards; ++i) {
//             bucket->tables[i] = g_hash_table_new_full(hash_func, key_equal_func,
//                                                       key_destroy_func, value_destroy_func);
//         }
//         bucket->flat = NULL;
//         bucket->shard_hash_func = hash_func;
// #else
//         fprintf(stderr, "wrong parameters for bucket, quitting ...\n");
//         exit(1);
// #endif
//     }
//     return bucket;
// }

GHashBucket* g_bucket_new_128(
    GDestroyNotify key_destroy_func,
    GDestroyNotify value_destroy_func,
    gsize          shards)
{
    GHashBucket *bucket = (GHashBucket*)g_malloc0(sizeof(GHashBucket));
    bucket->nshards = shards ? shards : 1;
    bucket->number_of_elements = 0;
    bucket->key_destroy_func = key_destroy_func;
    bucket->value_destroy_func = value_destroy_func;
    bucket->locks = (GMutex*)g_malloc0_n(bucket->nshards, sizeof(GMutex));
    for (gsize i = 0; i < bucket->nshards; ++i) g_mutex_init(&bucket->locks[i]);

    bucket->kind = BUCKET_FLAT128;
    bucket->flat = (FlatSet*)g_malloc0_n(bucket->nshards, sizeof(FlatSet));
    bucket->tables = NULL;
    bucket->shard_hash_func = NULL;

    return bucket;
}

void g_bucket_destroy(GHashBucket *b)
{
    if (!b) return;
#ifdef __USE_GLIB
    if (b->kind == BUCKET_GLIB) {
        for (gsize i = 0; i < b->nshards; ++i) g_hash_table_destroy(b->tables[i]);
        g_free(b->tables);
    } else
#endif
    {
        for (gsize i = 0; i < b->nshards; ++i) flat_free(&b->flat[i]);
        g_free(b->flat);
    }
    for (gsize i = 0; i < b->nshards; ++i) g_mutex_clear(&b->locks[i]);
    g_free(b->locks);
    g_free(b);
}

gsize g_bucket_size(GHashBucket *b) {
    return b ? ATOMIC_GET(b->number_of_elements) : 0;
}

gpointer g_bucket_lookup(GHashBucket* b, const key128_t* k)
{
    if (!b) return NULL;
#ifdef __USE_GLIB
    if (b->kind == BUCKET_GLIB) {
        size_t idx = shard_index_glib(b, key);
        G_MUTEX_LOCK(&b->locks[idx]);
        gpointer v = g_hash_table_lookup(b->tables[idx], key);
        G_MUTEX_UNLOCK(&b->locks[idx]);
        return v;
    } else
#endif
    {
        // const key128_t *k = (const key128_t*)key;
        size_t idx = shard_index_flat(b, k);
        G_MUTEX_LOCK(&b->locks[idx]);
        gboolean ok = flat_lookup(&b->flat[idx], k);
        G_MUTEX_UNLOCK(&b->locks[idx]);
        return ok ? GINT_TO_POINTER(TRUE) : NULL;
    }
}

gboolean g_bucket_remove(GHashBucket* b, const key128_t* k)
{
    if (!b) return FALSE;
    gboolean res = FALSE;
#ifdef __USE_GLIB
    if (b->kind == BUCKET_GLIB) {
        size_t idx = shard_index_glib(b, key);
        G_MUTEX_LOCK(&b->locks[idx]);
        if (g_hash_table_remove(b->tables[idx], key)) { b->number_of_elements--; res = TRUE; }
        G_MUTEX_UNLOCK(&b->locks[idx]);
    } else
#endif
    {
        // const key128_t *k = (const key128_t*)key;
        size_t idx = shard_index_flat(b, k);
        G_MUTEX_LOCK(&b->locks[idx]);
        res = flat_remove(&b->flat[idx], k);
        G_MUTEX_UNLOCK(&b->locks[idx]);
        if (res) {
            ATOMIC_DEC(b->number_of_elements);
        }
    }
    return res;
}

gboolean g_bucket_insert(GHashBucket *b, const key128_t* kptr, gpointer value)
{
    (void)value; if (!b) return FALSE;
 #ifdef __USE_GLIB   
    if (b->kind == BUCKET_GLIB) {
        size_t idx = shard_index_glib(b, key);
        G_MUTEX_LOCK(&b->locks[idx]);
        gboolean ins = g_hash_table_insert(b->tables[idx], key, GINT_TO_POINTER(TRUE));
        if (ins) b->number_of_elements++;
        G_MUTEX_UNLOCK(&b->locks[idx]);
        return ins;
    } else
#endif
    {
        // key128_t *kptr = (key128_t*)key;
        size_t idx = shard_index_flat(b, kptr);
        G_MUTEX_LOCK(&b->locks[idx]);
        //if (b->flat[idx].cap == 0) flat_init(&b->flat[idx], 1024);
        gboolean ins = flat_insert(&b->flat[idx], kptr);
        G_MUTEX_UNLOCK(&b->locks[idx]);
        /* We copied bytes; respect caller ownership per contract. */
        // if (b->key_destroy_func) b->key_destroy_func(key);
        if (ins) {
            ATOMIC_INC(b->number_of_elements);
        }
        return ins;
    }
}

gboolean g_bucket_insert_copy128(GHashBucket *b, const key128_t *key)
{
    if (!b) return FALSE;
#ifdef __USE_GLIB
    if (b->kind == BUCKET_GLIB) {
        key128_t *dup = (key128_t*)g_memdup2(key, sizeof(key128_t));
        size_t idx = shard_index_glib(b, dup);
        G_MUTEX_LOCK(&b->locks[idx]);
        gboolean ins = g_hash_table_insert(b->tables[idx], dup, GINT_TO_POINTER(TRUE));
        if (ins) b->number_of_elements++;
        else g_free(dup);
        G_MUTEX_UNLOCK(&b->locks[idx]);
        return ins;
    } else
#endif
    {
        size_t idx = shard_index_flat(b, key);
        G_MUTEX_LOCK(&b->locks[idx]);
        if (b->flat[idx].cap == 0) flat_init(&b->flat[idx], 1024);
        gboolean ins = flat_insert(&b->flat[idx], key);
        G_MUTEX_UNLOCK(&b->locks[idx]);
        if (ins) ATOMIC_INC(b->number_of_elements);
        return ins;
    }
}

void g_bucket_foreach128(GHashBucket *b, GKey128ForeachFunc func, gpointer user_data)
{
    if (!b || !func) return;
#ifdef __USE_GLIB
    if (b->kind != BUCKET_FLAT128) {
        /* Optional: warn; for GLIB we could iterate and cast keys if they are key128_t*,
         * but local sets in our use are always FLAT128. */
        return;
    }
#endif
    for (gsize s = 0; s < b->nshards; ++s) {
        G_MUTEX_LOCK(&b->locks[s]);
        FlatSet *fs = &b->flat[s];
        for (size_t i = 0; i < fs->cap; ++i) {
            if (fs->ctrl[i] == CTRL_FULL) func(&fs->keys[i], user_data);
        }
        G_MUTEX_UNLOCK(&b->locks[s]);
    }
}

gboolean g_bucket_is_flat128(GHashBucket* b) { return b && b->kind == BUCKET_FLAT128; }

// --- New: reserve total capacity across shards.
// NOTE: call this before concurrent use. No locking inside.
void g_bucket_reserve(GHashBucket *b, uint64_t total_expected)
{
    if (!b || b->kind != BUCKET_FLAT128 || total_expected == 0) return;

    // const double MAX_TRUE = 0.80; // docelowy load po rezerwacji
    const double SKEW     = 1.15; // zapas na nierówną dystrybucję; dostosuj w razie potrzeby

    // Prosty podział przez liczbę shardów + zapas na skew
    double base = (double) total_expected / (double) b->nshards;
    size_t per_shard_expected = (size_t)(base * SKEW);

    for (gsize i = 0; i < b->nshards; ++i) {
        flat_reserve(&b->flat[i], per_shard_expected, MAX_TRUE);
    }
}
