#include "common.h"
#include "bott.h"
#include "dag.h"

int main(void)
{
    const int dim = 9;
    vec_t *cache = NULL;
    size_t cache_size;

    printf("=== [bott-upper_triangular] testing dimension %d ===\n", dim);
    populate_cache(&cache, &cache_size, dim);

    state_t max_state = get_max_state(dim);

    vec_t *mat = init(dim);
    for (state_t s = 0; s <= max_state; ++s) {
        matrix_by_state(mat, cache, s, dim);
        if (! is_upper_triangular(mat, dim) ) {
            printf("    state %lu is not upper triangular:\n", s);
            print_mat(mat, dim);
            exit(1);
        }
    }
    free(mat);
    printf("=== [bott-upper_triangular] all tests passed ===\n");
    return 0;
}