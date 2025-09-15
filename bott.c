#include "bott.h"

void populate_cache(vec_t **cache, size_t *size, const ind_t dim)
{
    vec_t i, j, o, m;

    *size = 1<<(dim-2);
    if (*size<4) { *size = 4; }
    *cache = (vec_t*)malloc(VEC_T_SIZE*(*size));
    (*cache)[0] = 0;
    for (i=1, m=1; i<=dim-2; i++) {
        o = 3 << (dim-1-i);
        for (j=0; j < (1<<(i-1)); j++) {
            (*cache)[m++] = o ^ (*cache)[j];
        }
    }
}

vec_t *init(const ind_t dim)
{
    static vec_t *mat;
    mat = (vec_t*)malloc(VEC_T_SIZE * dim);
    memset(mat, 0, VEC_T_SIZE * dim);
    return mat;
}

void set(vec_t *mat, const vec_t *cache, const state_t state, const ind_t dim)
{
    vec_t c, j, i;
    state_t mask;
    for (i=1, c=1, j=dim-3; ; c+=i++, j--) {
        mask = ((state_t)1<<i) - 1;
        mat[j] = cache[(state>>(c-1))&mask];
        if (j==0) {
            break;
        }
    }
}

#ifdef DEBUG
void print(const vec_t *mat, const ind_t dim)
{
    ind_t i, j;
    static const char c[] = { '.', '1' };

    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            printf("%c", c[C(mat[i],dim,j)]);
        }
        printf("\n");
    }
}
#endif

/*
static inline vec_t scalar_product(const vec_t v1, const vec_t v2)
{
#if VEC_T_SIZE == 64
    return __builtin_parityl(v1&v2);
#else
    return __builtin_parity(v1&v2);
#endif
}
*/

static inline int equal_cols(const vec_t *mat, const ind_t dim, const ind_t i, const ind_t j)
{
    vec_t mask = (1<<i) ^ (1<<j);
    ind_t r;
    for (r=0; r<dim; r++) {
        if (scalar_product(mat[r],mask)) {
            return 0;
        }
    }
    /* r = i; */
    // for (; r<j; r++) {
    //     if (mat[r]&mask) {
    //         return 0;
    //     }
    // }
    return 1;
}

size_t is_spinc(const vec_t *mat, const ind_t dim)
{
    ind_t i, j, e, z;
    vec_t aij;
    for (j=2; j<dim-2; j++) {
        //e = 1;
        //z = 1;
        for (i=0, e=1, z=1; i<j; i++) {
            /*
            if (equal_cols(mat, dim, i, j)) {
                aij = 0;
            } else {
                aij = scalar_product(mat[i], mat[j]);
            }
            */
            aij = equal_cols(mat, dim, i, j) ? 0 : scalar_product(mat[i], mat[j]);
            if (z && aij) {
                z = 0;
            }
            if (e && aij!=C(mat[i],dim,j)) {
                e = 0;
            }
            if (e==0 && z==0) {
                return 0;
            }
        }
    }
    return 1;
}

static inline int row_sum(vec_t r)
{
#if VEC_T_SIZE == 64
    return __builtin_popcountl(r);
#else
    return __builtin_popcount(r);
#endif
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
