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
    fprintf(stderr, "-s: starting dimension, between 3 and 11\n");
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

int calculate_spin = 0;

/*
 * cdim - current dimension
 * ddim - destination dimension
*/
size_t backtrack(vec_t *mat, const vec_t *cache, vec_t *tr, ind_t cdim, ind_t ddim, size_t *spinc, size_t *spin)
{
    vec_t max = 1<<(cdim-2);
    vec_t r;
    vec_t row = ddim-cdim;

    if (row==0) {
        /* 
         * final step:
        * 
        * no need to calculate matrix 
        * [ 0 000 ]
        * [ 0 mat ]
        * since spinc info from mat
        * is the same
        */
        mat[0]  = 0;
        *spinc += 1;
        if (calculate_spin) {
            *spin  += is_spin(mat, ddim);
            for (r = 1; r<max; r++) {
                mat[0] = cache[r];
                transpose(mat, tr, ddim, 1);
                if (is_spinc_tr(mat,ddim, tr)) {
                    *spinc += 1;
                    *spin  += is_spin(mat, ddim);
                }
            }
        } else {
            for (r = 1; r<max; r++) {
                mat[0] = cache[r];
                *spinc += is_spinc(mat, ddim);
            }
        }
    }
    else {
        /* recursive call */
        for (r=0; r<max; r++) {
            mat[row] = cache[r];
            transpose(&mat[row], &tr[row], cdim, 1);
            if (is_spinc_tr(&mat[row], cdim, &tr[row])) {
                backtrack(mat, cache, tr, cdim+1, ddim, spinc, spin);
            }
        }
    }
    return *spinc;
}

int main(int argc, char *argv[])
{
    int opt = 1, v = 0;
    long time;
    /* default values of dimension, num of threads */
    ind_t sdim = 4, dim = 6;
    vec_t *mat, *cache = NULL, *tr;
    /* state is a number which is used to generate RBM matrix */
    vec_t state, max_state, row;
    size_t spinc, spin, a, b;
    size_t cache_size;
    while ((opt = getopt(argc, argv, "vhj:d:s:a")) != -1) {
        switch (opt) {
        case 'a':
            calculate_spin = 1;
            break;
        case 'v':
            v = 1;
            break;
        case 'j':
            omp_set_num_threads(atoi(optarg));
            break;
        case 'd':
            dim = atoi(optarg);
            break;
        case 's':
            sdim = atoi(optarg);
            if (sdim < 3 || sdim > 11) {
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

    tic();
    row = dim-sdim;
    spinc = 0;
    spin  = 0;
    #pragma omp parallel default(none) private(mat, a, b, tr) shared(cache, max_state, dim, sdim, row) reduction(+:spinc,spin)
    {
        a = 0;
        b = 0;    
        mat = init(dim);
        tr  = init(dim);
        #pragma omp for
        for (state=0; state<=max_state; state++) {
            set(&mat[row], cache, state, sdim);
            transpose(&mat[row], &tr[row], sdim, sdim);
            if (is_spinc_tr(&mat[row], sdim, &tr[row])) {
                backtrack(mat, cache, tr, sdim+1, dim, &a, &b);
            }
        }
        free(mat);
        spinc+= a;
        spin += b;
    }
    time = toc();

    free(cache);

    printf("spinc: %lu\n", spinc);
    if (calculate_spin) {
        printf("spin:  %lu\n", spin);
    }
    if (v) {
        printf("time:  %.03fs\n", (float)time/1000000000);
    }

    return 0;
}

