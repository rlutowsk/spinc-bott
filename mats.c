#include <omp.h>

#include "bott.h"
#include "dag.h"

void help(const char *name)
{
    fprintf(stderr, "Usage: %s [-d dimension] [-j njobs] [-v] [-h] [-p]\n", name);
}

int main(int argc, char *argv[])
{
    int opt = 1, v = 0;
    long time;
    /* default values of dimension, num of threads */
    ind_t dim = 6, j=16, c;
    vec_t *cache = NULL;
    /* state is a number which is used to generate RBM matrix */
    state_t state, max_state;
    size_t spinc, spin;
    size_t cache_size;
    int p = 0;
    while ((opt = getopt(argc, argv, "vhj:d:p")) != -1) {
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
        case 'p':
            p = 1;
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
    tic();
    #pragma omp parallel shared(cache, dim, max_state) reduction(+:spin,spinc)
    {
        vec_t *mat   = init(dim);
        size_t a = 0, b = 0;

        char d6_buf[128];
        size_t cnt = 0;
        size_t siz = strlen( matrix_to_d6(mat, dim, d6_buf) );
        char *buf = (char*)malloc( (siz+1)*10000+1 );
        if (buf==NULL) {
            fprintf(stderr, "Error in malloc\n");
            exit(1);
        }
        char *pnt = buf;
        
        #pragma omp for
        for (state=0; state<=max_state; state++) {
            matrix_by_state(mat, cache, state, dim);
            if (is_spinc(mat, dim)) {
                a++;
                b += is_spin(mat, dim);
                if (p) {
                    matrix_to_d6(mat, dim, d6_buf);
                    memcpy( pnt, d6_buf, siz);
                    pnt += siz;
                    *pnt++ = '\n';
                    if (++cnt == 10000) {
                        cnt = 0;
                        *pnt = '\0';
                        pnt = buf;
                        #pragma omp critical
                        {
                            printf("%s", buf);
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            *pnt = '\0';
            printf("%s", buf);
        }

        free(buf);
        free(mat);
        spinc += a;
        spin  += b;
    }
    time = toc();

    free(cache);

    fprintf(stderr, "spinc: %lu\n", spinc);
    fprintf(stderr, "spin:  %lu\n", spin);
    if (v) {
        fprintf(stderr, "time:  %.03fs\n", (float)time/1000000000);
    }

    return 0;
}

