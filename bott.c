#include <assert.h>

#include "bott.h"

void populate_cache(vec_t **cache, size_t *size, const int dim)
{
    *size = 1<<(dim-2);
    if (*size<4) { *size = 4; }
    *cache = (vec_t*)malloc(sizeof(vec_t)*(*size));
    // read from right to left:
    (*cache)[0] = 0;          // 000 - dim 2,3,4
    (*cache)[1] = 3<<(dim-2); // 011 - dim   3,4
    (*cache)[2] = 3<<(dim-3); // 110 - dim     4
    (*cache)[3] = 5<<(dim-3); // 101 - dim     4
    size_t idx = 4;
    for (int i = 3; i < dim-1; ++i) {
        // additional basis vector:
        vec_t new_vec = (vec_t)3 << (dim - i - 1);
        // add the newly created vector to all previously created vectors:
        for (int j = 0; j < (1 << (i - 1)); ++j, ++idx) {
            (*cache)[idx] = (*cache)[j] ^ new_vec;
        }
    }
    assert(idx == *size);
}

vec_t *init(const ind_t dim)
{
    vec_t *mat = (vec_t*)calloc(dim, sizeof(vec_t));
    return mat;
}

void print_mat(const vec_t *mat, const ind_t dim)
{
    static const char c[] = { '.', '1' };

    for (ind_t i=0; i<dim; i++) {
        for (ind_t j=0; j<dim; j++) {
            printf("%c", c[C(mat[i],j)]);
        }
        printf("\n");
    }
}

static INLINE bool equal_cols(const vec_t *mat, const ind_t dim, const ind_t i, const ind_t j)
{
    vec_t mask = (1<<i) ^ (1<<j);
    for (ind_t r=0; r<dim; r++) {
        if (scalar_product(mat[r],mask)) {
            return false;
        }
    }
    return true;
}

bool is_spinc(const vec_t *mat, const ind_t dim)
{
    if (!is_orientable(mat, dim)) {
        return false;
    }
    for (ind_t j=2; j<dim-2; j++) {
        for (ind_t i=0, e=1, z=1; i<j; i++) {
            vec_t aij = equal_cols(mat, dim, i, j) ? 0 : scalar_product(mat[i], mat[j]);
            if (z && aij) {
                z = 0;
            }
            if (e && aij!=C(mat[i],j)) {
                e = 0;
            }
            if (e==0 && z==0) {
                return false;
            }
        }
    }
    return true;
}

bool is_spin(const vec_t *mat, const ind_t dim)
{
    if (!is_orientable(mat, dim)) {
        return false;
    }
    for (ind_t j=1; j<dim-2; j++) {
        int rem = row_sum(mat[j]) % 4;
        if (rem == 2) {
            for (ind_t i=0; i<j; i++) {
                if (scalar_product(mat[i], mat[j]) != C(mat[i],j)) {
                    return false;
                }
            }
        } else {
            for (ind_t i=0; i<j; i++) {
                if (scalar_product(mat[i], mat[j])==1) {
                    return false;
                }
            }
        }
    }
    return true;
}

/*
 * isomorphism operations for Bott matrices
 * based on:
 * "Real Bott manifolds and acyclic digraphs" by Suyoung Choi, Mikiya Masuda, Sang-il Oum
 * https://arxiv.org/abs/1104.0072
 */

/* Op1 */
void swap_rows_and_cols(vec_t *src, vec_t *dst, ind_t dim, ind_t r1, ind_t r2)
{
    vec_t diff, mask;

    for (ind_t i=0; i<dim; i++) {
        diff = ((src[i] >> r1) ^ (src[i] >> r2)) & 1;
        mask = (diff << r1) | (diff << r2);
        dst[i] = src[i] ^ mask;
    }
    SWAP(vec_t, dst[r1], dst[r2]);
}

/* Op2 */
void conditional_add_col(const vec_t *src, vec_t *dst, ind_t dim, ind_t k)
{
    vec_t row_k = src[k];

    for (ind_t i=0; i<dim; i++) {
        vec_t cond = -(vec_t)(((src[i] >> k) & 1));
        dst[i] = src[i] ^ (row_k & cond);
    }
}

/* Op3 */
bool conditional_add_row(const vec_t *src, vec_t *dst, ind_t dim, ind_t l, ind_t m)
{
    if ( (l==m) || ! (equal_cols(src, dim, l, m)) ) {
        return false;
    }

    for (ind_t i=0; i<dim; i++) {
        dst[i] = src[i];
    }
    dst[m] = src[m] ^ src[l];

    return true;
}

/* return the numbers of ones in the matrix */
int matrix_weight(const vec_t *mat, const ind_t dim)
{
    int w = 0;
    for (ind_t i=0; i<dim; i++) {
        w += row_sum(mat[i]);
    }
    return w;
}
