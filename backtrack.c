#include "bott.h"

/* measurements */
#include <time.h>
struct timespec tv;
void tic(void)
{
    clock_gettime(CLOCK_REALTIME, &tv);
}
long toc(void)
{
    static struct timespec st;
    clock_gettime(CLOCK_REALTIME, &st);
    return (tv.tv_sec==st.tv_sec) ? st.tv_nsec - tv.tv_nsec : 1000000000*(st.tv_sec-tv.tv_sec)+st.tv_nsec - tv.tv_nsec;
}

void help(const char *name)
{
    fprintf(stderr, "Usage: %s [-d dimension] [-j njobs] [-s start_dim] [-v] [-h]\n", name);
    fprintf(stderr, "-d: dimension to calculate\n");
    fprintf(stderr, "-j: number of threads\n");
    fprintf(stderr, "-s: starting dimension, between 3 and 9\n");
    fprintf(stderr, "-v: be verbose\n");
}

inline vec_t get_max_state(ind_t dim)
{
    ind_t c;
    vec_t max_state;
    for (c=0, max_state=1; c<dim-1; c++) {
        max_state <<= (dim-c-2);
    }
    return max_state - 1;
}

struct result {
    size_t spinc;
    size_t spin;
};

/*
 * cdim - current dimension
 * ddim - destination dimension
*/
size_t backtrack(vec_t *mat, const vec_t *cache, ind_t cdim, ind_t ddim)
{
    size_t spinc = 0;
    vec_t max = 1<<(cdim-2);
    vec_t r;
    vec_t row;

    if (cdim==ddim) {
        /* final step */
        for (r = 0; r<max; r++) {
            mat[0] = cache[r];
            if (is_spinc(mat, ddim)) {
                spinc++;
            }
        }
    }
    else {
        /* recursive call */
        row = ddim-cdim;
        for (r=0; r<max; r++) {
            mat[row] = cache[r];
            if (is_spinc(&mat[row], cdim)) {
                spinc += backtrack(mat, cache, cdim+1, ddim);
            }
        }
    }
    return spinc;
}

int main(int argc, char *argv[])
{
    int opt = 1, v = 0;
    long time;
    /* default values of dimension, num of threads */
    ind_t sdim = 4, dim = 6, j=16;
    vec_t *mat, *cache = NULL;
    /* state is a number which is used to generate RBM matrix */
    vec_t state, max_state;
    size_t spinc, spin;
    size_t cache_size;
    while ((opt = getopt(argc, argv, "vhj:d:s:")) != -1) {
        switch (opt) {
        case 'v':
            v = 1;
            break;
        case 'j':
            j = atoi(optarg);
            omp_set_num_threads(j);
            break;
        case 'd':
            dim = atoi(optarg);
            break;
        case 's':
            sdim = atoi(optarg);
            if (sdim < 3 || sdim > 10) {
                help(argv[0]);
                exit(EXIT_FAILURE);
            }
            break;
        case 'h':
            help(argv[0]);
            exit(EXIT_SUCCESS);
        default: /* '?' */
            help(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (sdim >= dim) {
        fprintf(stderr, "Starting dimension must be less than the destination dimension.\n");
        exit(EXIT_FAILURE);
    }

    populate_cache(&cache, &cache_size, dim);

    max_state = get_max_state(sdim);

    spinc = 0;
    tic();
    #pragma omp parallel private (mat) shared(cache, dim, max_state)
    {
        spinc = 0;
        spin  = 0;
        mat   = init(dim);
        #pragma omp for reduction(+:spinc,spin) //schedule(static)
        for (state=0; state<=max_state; state++) {
            set(mat, cache, state, dim);
            if (is_spinc(mat, dim)) {
                spinc += backtrack(mat, cache, sdim+1, dim);
            }
        }
        free(mat);
    }
    time = toc();

    free(cache);

    printf("spinc: %lu\n", spinc);
    printf("spin:  %lu\n", spin);
    if (v) {
        printf("time:  %.03fs\n", (float)time/1000000000);
    }

    return 0;
}

