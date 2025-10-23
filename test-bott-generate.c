#include "common.h"
#include "bott.h"
#include "dag.h"
#include "bucket.h"
#include "adjpack11.h"

int main(void)
{
    const int dim = 8;
    vec_t *cache = NULL;
    size_t cache_size;
    populate_cache(&cache, &cache_size, dim);

    state_t max_state = get_max_state(dim);

    GHashBucket *code_set = g_bucket_new_128(NULL, NULL, 1023);

    vec_t *mat = init(dim);
    key128_t k;
    printf("=== [bott-generate] testing dimension %d for duplicate codes ===\n", dim);
    for (state_t s = 0; s <= max_state; ++s) {
        if ( g_bucket_size(code_set) != s ) {
            fprintf(stderr, "    code duplicates error at state %lu\n", s-1);
            exit(1);
        }
        matrix_by_state(mat, cache, s, dim);
        adjpack_from_matrix(mat, dim, &k);
        g_bucket_insert_copy128(code_set, &k);
    }
    printf("    Found %lu unique codes for dimension %d, as expected.\n", g_bucket_size(code_set), dim);
    printf("=== [bott-generate] all tests passed ===\n");
    free(mat);

    g_bucket_destroy(code_set);

    free(cache);

    return 0;
}
