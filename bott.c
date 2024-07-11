#include "bott.h"

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

/* set transpose in s <= row < e */
void transpose(const vec_t *src, vec_t *dst, const ind_t dim, ind_t e)
{
    ind_t i,j;
    /* first column is always zero */
    dst[0] = 0;
    for (j=1; j<e; j++) {
        dst[j] = 0;
        for (i=0; i<j; i++) {
            dst[j] |= (C(src[i], dim, j))<<(dim-i-1);
        }
    }
    for (j=e; j<dim; j++) {
        dst[j] &= (1<<(dim-e))-1;
        for (i=0; i<e; i++) {
            dst[j] |= (C(src[i], dim, j))<<(dim-i-1);
        }
    }
}

static inline int equal_cols(const vec_t *mat, const ind_t dim, const ind_t i, const ind_t j)
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

static inline vec_t scalar_product(const vec_t v1, const vec_t v2)
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

static inline int equal_cols_tr(const vec_t *tr, ind_t i, ind_t j)
{
    return tr[i]==tr[j];
}

size_t is_spinc_tr(const vec_t *mat, const ind_t dim, const vec_t *tr)
{
    ind_t i, j, e, z;
    vec_t aij;
    for (j=2; j<dim-2; j++) {
        e = 1;
        z = 1;
        for (i=0; i<j; i++) {
            /*
            if (equal_cols(mat,dim,i,j)!=(tr[i]==tr[j])) {
                return 0;
            }
            */
            if (equal_cols_tr(tr,i,j)) {
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

static inline int row_sum(vec_t r)
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

