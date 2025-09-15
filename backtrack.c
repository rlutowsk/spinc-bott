#include "bott.h"
#include "common.h"

void help(const char *name)
{
    fprintf(stderr, 
        "Usage: %s [-j njobs] [-s start_dim] [-d dimension] [-v] [-h]\n" \
        "-j: number of threads\n" \
        "-d: target dimension to calculate\n" \
        "-s: starting dimension, between 3 and 11, less than target dimension\n" \
        "-a: calculate spin structures also\n" \
        "-v: be verbose\n",
        name );
}

static inline state_t get_max_state(ind_t dim)
{
    ind_t c;
    state_t max_state;
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
size_t backtrack(vec_t *mat, vec_t **cache, ind_t cdim, ind_t ddim, size_t *spinc, size_t *spin);

int main(int argc, char *argv[])
{
    int opt = 1, bs=-1, no_output=0;
    state_t test = 0;
#ifdef RBM_SHOW_PROGRESS
    int progress=0;
#endif
    /* default values of dimension, num of threads */
    ind_t sdim = 0, dim = 6;
    vec_t *mat, *cache[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    /* state is a number which is used to generate RBM matrix */
    state_t state, max_state;
    state_t loop_start=1, loop_stop=1;
    vec_t row;
    size_t spinc, spin, a, b;
    size_t cache_size;

    tic();

    while ((opt = getopt(argc, argv, "vhj:d:s:ac:npt:")) != -1) {
        switch (opt) {
        case 't':
            test = atol(optarg);
            break;
		case 'p':
        #ifdef RBM_SHOW_PROGRESS
			progress = 1;
        #else
            printlog(0, "compiled without progress support, flag '-p' not effective\n");
        #endif
			break;
        case 'n':
            no_output = 1;
            break;
        case 'c':
            bs = atoi(optarg);
            break;
        case 'a':
            calculate_spin = 1;
            break;
        case 'v':
            verbosity_level++;
            break;
        case 'j':
            omp_set_num_threads(atoi(optarg));
            break;
        case 'd':
            dim = (ind_t)atoi(optarg);
            break;
        case 's':
            sdim = (ind_t)atoi(optarg);
            break;
        case 'h':
            help(argv[0]);
            exit(EXIT_SUCCESS);
        default: /* '?' */
            help(argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    if (sdim==0) {
        sdim = (dim>11)? 11 : dim-1;
    }
    if (sdim < 3 || sdim > 11 || sdim>=dim) {
        fprintf(stderr, "Starting dimension out of range, see help.\n");
        exit(EXIT_FAILURE);
    }

    for (int i=sdim; i<12; i++) {
        vec_t *tmp;
        populate_cache(&tmp, &cache_size, i);
        cache[i] = tmp;
    }

    max_state = get_max_state(sdim);

    row = dim-sdim;
    spinc = 0;
    spin  = 0;

    /* try to guess proper chunk size */
    if (bs==-1) {
        bs = 7*sdim-49;
        bs = (bs<0)?0:bs;
    }

    if (test>0) {
        loop_start = (max_state+1)/2;
        loop_stop  = loop_start + (1<<test) - 1;
    } else {
        loop_start = 0;
        loop_stop  = max_state;
    }
    printlog(2, "start: %lu, end: %lu\n", loop_start, loop_stop);
    omp_set_schedule(omp_sched_dynamic, 1<<bs);
    printlog(2, "openmp chunk size set to 2^%lu\n", bs);

    #pragma omp parallel private(mat, a, b) shared(bs, stdout, cache, max_state, dim, sdim, row, calculate_spin, loop_start, loop_stop) reduction(+:spinc,spin)
    {
    #ifdef RBM_SHOW_PROGRESS
		vec_t p = 0, sp = (bs>0)?1<<(bs-1):1, total=max_state+1;
		int show_progress = progress && (omp_get_thread_num()==0);
        time_t x = toc(), y;
        double percent;
    #endif
        
		a = 0;
        b = 0;    
        mat = init(dim);
        #pragma omp for nowait schedule(runtime)
        for (state=loop_start; state<=loop_stop; state++) {
        #ifdef RBM_SHOW_PROGRESS    
			if (show_progress)
			{
				if (p==0 && (y=toc())>x+1000000000) {
                    x = y;
                    percent = 100.0*state/total;
					fprintf(stdout, "\33[2K\r%8.4f%% (ETT: %.1fs)\r", percent, 1.0*x/(percent*10000000));
					fflush(stdout);
				}
				if (++p == sp) {
					p = 0;
				}
			}
        #endif
            set(&mat[row], cache[sdim], state, sdim);
            if (is_spinc(&mat[row], sdim)) {
                backtrack(mat, cache, sdim+1, dim, &a, &b);
            }
        }
        free(mat);
        #pragma omp critical 
        {
            if (calculate_spin) {
                printlog(2, "thread %2d finished: %lu/%lu spin/spinc manifolds\n", omp_get_thread_num(), b, a);
            } else {
                printlog(2, "thread %2d finished: %lu spinc manifolds\n", omp_get_thread_num(), a);
            }
        }
        spinc+= a;
        spin += b;
    }

    if (calculate_spin) {
        printlog(1, "calculations finished: %lu/%lu spin/spinc manifolds in total\n", spin, spinc);
    } else {
        printlog(1, "calculations finished: %lu spinc manifolds in total\n", spinc);
    }

    if (!no_output) {
        printf("\33[2K\r");
        if (calculate_spin) {
            printf("%lu/%lu\n", spin, spinc);
        } else {
            printf("%lu\n", spinc);
        }
    }

    for(int i=sdim; i<12; i++) {
        free(cache[i]);
    }

    return 0;
}

static inline void increase_dimension(vec_t *mat, ind_t dim)
{
    for (ind_t i=1; i<dim; i++) {
        mat[i] <<= 1;
    }
}

static inline void decrease_dimension(vec_t *mat, ind_t dim)
{
    for (ind_t i=1; i<dim; i++) {
        mat[i] >>= 1;
    }
}

size_t backtrack(vec_t *mat, vec_t **cache, ind_t cdim, ind_t ddim, size_t *spinc, size_t *spin)
{
    vec_t max = 1<<(cdim-2);
    vec_t r;
    vec_t row = ddim-cdim;

    increase_dimension(mat, ddim);

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
                mat[0] = cache[ddim][r];
                if (is_spinc(mat,ddim)) {
                    *spinc += 1;
                    *spin  += is_spin(mat, ddim);
                }
            }
        } else {
            for (r = 1; r<max; r++) {
                mat[0] = cache[ddim][r];
                *spinc += is_spinc(mat, ddim);
            }
        }
    }
    else {
        /* recursive call */
        for (r=0; r<max; r++) {
            mat[row] = cache[cdim][r];
            if (is_spinc(&mat[row], cdim)) {
                backtrack(mat, cache, cdim+1, ddim, spinc, spin);
                /* 
                 * after every call we have to decrease dimension of all rows
                 * since next call will increase it again
                 */
                decrease_dimension(mat, ddim);
            }
        }
    }
    return *spinc;
}
