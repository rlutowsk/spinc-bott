#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

/* ensure we are in 64-bit environment */
#include <bits/wordsize.h>
#if __WORDSIZE != 64
# error "compile only in 64-bit mode"
#endif

#define SWAP(T, a, b) do { T SWAP_TMP = a; a = b; b = SWAP_TMP; } while (0)
#define C(row,dim,j) ((row>>(dim-j-1))&1)

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

typedef unsigned char ind_t;
typedef unsigned long vec_t;

void help(const char *name)
{
    fprintf(stderr, "Usage: %s [-d dimension] [-j njobs] [-v] [-h]\n", name);
}

/* cache holds possible values for bit-sequences of orientable RBMs */
void populate_cache(vec_t **cache, size_t *size, const ind_t dim);

/* allocate memory for RBM matrix */
vec_t *init(const ind_t dim);

/* set values of rows of RBM matrix dpending on state */
void set(vec_t *mat, const vec_t *cache, const vec_t state, const ind_t dim);

/* main workers of this app */
size_t is_spinc(const vec_t *mat, const ind_t dim);
size_t is_spin(const vec_t *mat, const ind_t dim);

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

void populate_cache(vec_t **cache, size_t *size, const ind_t dim)
{
    vec_t i, j, o, m;

    *size = 1<<(dim-2);
    if (*size<4) { *size = 4; }
    *cache = (vec_t*)malloc(__WORDSIZE*(*size));
    (*cache)[0] = 0;
    (*cache)[1] = 3;
    (*cache)[2] = 5;
    (*cache)[3] = 6;
    for (i=4, m=6; i<*size; ) {
        o = i>>1;
        m <<= 1;
        for (j=i; j<i+o; j++) {
            (*cache)[j] = m^(*cache)[j-o];
        }
        for (i<<=1, o *= 3; j<i; j++) {
            (*cache)[j] = m^(*cache)[j-o];
        }
    }

}

vec_t *init(const ind_t dim)
{
    static vec_t *mat;
    mat = (vec_t*)malloc(__WORDSIZE*dim);
    bzero(mat, __WORDSIZE*dim);
    return mat;
}

void set(vec_t *mat, const vec_t *cache, const vec_t state, const ind_t dim)
{
    vec_t c, j, i, mask;
    for (i=1, c=1, j=dim-3; ; c+=i++, j--) {
        mask = (1<<i) - 1;
        mat[j] = cache[(state>>(c-1))&mask];
        if (j==0) {
            break;
        }
    }
}

inline int equal_cols(const vec_t *mat, const ind_t dim, const ind_t i, const ind_t j)
{
    vec_t mask = (1<<(dim-i-1)) ^ (1<<(dim-j-1));
    ind_t r;
    for (r=0; r<j; r++) {
        if (__builtin_parityl((mat[r])&mask)) {
            return 0;
        }
    }
    return 1;
}

inline vec_t scalar_product(const vec_t v1, const vec_t v2)
{
    return __builtin_parityl(v1&v2);
}

size_t is_spinc(const vec_t *mat, const ind_t dim)
{
    ind_t i, j, e, z;
    vec_t aij;
    for (j=2; j<dim-2; j++) {
        e = 1;
        z = 1;
        for (i=0; i<j; i++) {

            if (equal_cols(mat, dim, i, j)) {
                aij = 0;
            } else {
                aij = scalar_product(mat[i], mat[j]);
            }
            if (z && aij) {
                z = 0;
            }
            if (e && aij!=C(mat[i],dim,j)) {
                e = 0;
            }
            if (!e && !z) {
                return 0;
            }
        }
    }
    return 1;
}

inline int row_sum(vec_t r)
{
    return __builtin_popcountl(r);
}

size_t is_spin(const vec_t *mat, const ind_t dim)
{
    ind_t i,j;
    for (j=1; j<dim-2; j++) {
        if (row_sum(j)%4 == 2) {
            for (i=0; i<j; i++) {
                if (scalar_product(mat[i], mat[j]) != C(mat[i],dim,j)) {
                    return 0;
                }
            }
        } else {
            for (i=0; i<j; i++) {
                if (scalar_product(mat[i], mat[j])==1) {
                    return 0;
                }
            }
        }
    }
    return 1;
}

