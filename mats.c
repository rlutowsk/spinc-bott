#include "bott.h"
#include "common.h"

void help(const char *name)
{
    fprintf(stderr, "Usage: %s [-d dimension] [-j njobs] [-v] [-h]\n", name);
}

int main(int argc, char *argv[])
{
    int opt = 1, v = 0;
    long time;
    /* default values of dimension, num of threads */
    ind_t dim = 6, j=16, c;
    vec_t *cache = NULL;
    mat_t mat = 0;
    /* state is a number which is used to generate RBM matrix */
    state_t state, max_state;
    size_t spinc, spin, a, b;
    size_t cache_size;
    while ((opt = getopt(argc, argv, "vhj:d:")) != -1) {
        switch (opt) {
        case 'v':
            v = 1;
            break;
        case 'j':
            j = (ind_t)atoi(optarg);
            omp_set_num_threads(j);
            break;
        case 'd':
            dim = (ind_t)atoi(optarg);
            break;

        case 'h':
            help(argv[0]);
            exit(EXIT_SUCCESS);
        default: /* '?' */
            help(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    for (c=0, max_state=1; c<dim-1; c++) {
        max_state <<= (dim-c-2);
    }
    max_state -= 1;

    populate_cache(&cache, &cache_size, dim);

    spinc = 0;
    spin  = 0;
    new_init(dim);
    tic();
    #pragma omp parallel private (mat, a, b) shared(cache, dim, max_state) reduction(+:spin,spinc)
    {
        // old_mat   = init(dim);
        a = 0;
        b = 0;
        #pragma omp for
        for (state=0; state<=max_state; state++) {
            new_set(&mat, cache, state, dim);
            if (new_is_spinc(mat, dim)) {
                a++;
                //b += is_spin(mat, dim);
            }
        }
        // free(old_mat);
        spinc += a;
        spin  += b;
    }
    time = toc();

    free(cache);

    printf("spinc: %lu\n", spinc);
    //printf("spin:  %lu\n", spin);
    if (v) {
        printf("time:  %.03fs\n", (float)time/1000000000);
    }

    return 0;
}

