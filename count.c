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
    fprintf(stderr, "Usage: %s [-d dimension] [-j njobs] [-v] [-h]\n", name);
}

int main(int argc, char *argv[])
{
    int opt = 1, v = 0;
    long time;
    /* default values of dimension, num of threads */
    ind_t dim = 6, j=16, c;
    vec_t *mat, *cache = NULL;
    /* state is a number which is used to generate RBM matrix */
    vec_t state, max_state;
    size_t spinc, spin;
    size_t cache_size;
    while ((opt = getopt(argc, argv, "vhj:d:")) != -1) {
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
    tic();
    #pragma omp parallel private (mat) shared(cache, dim, max_state)
    {
        spinc = 0;
        spin  = 0;
        mat   = init(dim);
        #pragma omp for reduction(+:spinc,spin) schedule(static)
        for (state=0; state<=max_state; state++) {
            set(mat, cache, state, dim);
            if (is_spinc(mat, dim)) {
                spinc++;
                spin += is_spin(mat, dim);
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
        printf("max:   %lu\n", max_state+1);
        printf("rate:  %.03f\n", (float)spinc/(max_state+1)*100);
    }

    return 0;
}

