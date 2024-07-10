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

static inline vec_t get_max_state(ind_t dim)
{
    ind_t c;
    vec_t max_state;
    for (c=0, max_state=1; c<dim-1; c++) {
        max_state <<= (dim-c-2);
    }
    return max_state - 1;
}

/*
 * cdim - current dimension
 * ddim - destination dimension
*/
size_t backtrack(vec_t *mat, const vec_t *cache, ind_t cdim, ind_t ddim, size_t *spinc, size_t *spin)
{
    vec_t max = 1<<(cdim-2);
    vec_t r;
    vec_t row;

    if (cdim==ddim) {
        /* final step */
        for (r = 0; r<max; r++) {
            mat[0] = cache[r];
            if (is_spinc(mat, ddim)) {
                *spinc += 1;
                *spin += is_spin(mat, ddim);
            }
        }
    }
    else {
        /* recursive call */
        row = ddim-cdim;
        for (r=0; r<max; r++) {
            mat[row] = cache[r];
            if (is_spinc(&mat[row], cdim)) {
                backtrack(mat, cache, cdim+1, ddim, spinc, spin);
            }
        }
    }
    return *spinc;
}

size_t backtrack1(vec_t *mat, const vec_t *cache, ind_t cdim, ind_t pdim, ind_t ddim, size_t *spinc, size_t *spin)
{
    vec_t max = 1<<(cdim-2);
    vec_t r, i;
    vec_t row;
    vec_t *tmat;
    size_t _spin, _spinc, a, b;

    if (cdim > pdim) {
        return 0;
    }
    _spinc = 0;
    _spin  = 0;
    row = ddim-cdim;
    if (cdim == pdim) {
        #pragma omp parallel default(none) private(tmat, a, b, i) shared(ddim, row, max, cache, cdim, mat) reduction(+:_spinc,_spin)
        {
            a = 0;
            b = 0;
            tmat = init(ddim);
            for (i=0; i<row; i++) {
                tmat[i] = 0;
            }
            for (i=row+1; i<ddim; i++) {
                tmat[i] = mat[i];
            }
            #pragma omp for schedule(static)
            for (r=0; r<max; r++) {
                tmat[row] = cache[r];
                if (is_spinc(&tmat[row], cdim)) {
                    backtrack(tmat, cache, cdim+1, ddim, &a, &b);
                }
            }
            free(tmat);
            _spinc += a;
            _spin  += b;
        }
        *spinc += _spinc;
        *spin  += _spin;

        return *spinc;
    }
    /* recursive call */
    for (r=0; r<max; r++) {
        mat[row] = cache[r];
        if (is_spinc(&mat[row], cdim)) {
            backtrack1(mat, cache, cdim+1, pdim, ddim, spinc, spin);
        }
    }
    
    return *spinc;
}

int main(int argc, char *argv[])
{
    int opt = 1, v = 0;
    long time;
    /* default values of dimension, num of threads */
    ind_t sdim = 4, dim = 6, pdim = 5, j=16;
    vec_t *mat, *cache = NULL;
    /* state is a number which is used to generate RBM matrix */
    vec_t state, max_state, row;
    size_t spinc, spin, bits;
    size_t cache_size;
    while ((opt = getopt(argc, argv, "vhj:d:s:p:")) != -1) {
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
        case 'p':
            pdim = atoi(optarg);
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
    if (pdim<=sdim || pdim>=dim) {
        fprintf(stderr, "Parallel spit dimension must be in bounds [%d, %d).\n", sdim, dim);
        exit(EXIT_FAILURE);
    }

    populate_cache(&cache, &cache_size, dim);

    max_state = get_max_state(sdim);

    tic();
    row = dim-sdim;
    bits= __WORDSIZE*row;
    spinc = 0;
    spin  = 0;
    mat   = init(dim);
    for (state=0; state<=max_state; state++) {
        memset(mat, 0, bits);
        set(&mat[row], cache, state, sdim);
        if (is_spinc(&mat[row], sdim)) {
            backtrack1(mat, cache, sdim+1, pdim, dim, &spinc, &spin);
        }
    }
    free(mat);

    time = toc();

    free(cache);

    printf("spinc: %lu\n", spinc);
    printf("spin:  %lu\n", spin);
    if (v) {
        printf("time:  %.03fs\n", (float)time/1000000000);
    }

    return 0;
}

