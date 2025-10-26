#ifndef BOTT_H
#define BOTT_H

#include "common.h"

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

/* ensure we are in 64-bit environment */
#include <bits/wordsize.h>
#if __WORDSIZE != 64
    #error "compile only in 64-bit mode"
#endif

#ifndef SWAP
    #define SWAP(T, a, b) do { T SWAP_TMP = a; a = b; b = SWAP_TMP; } while (0)
#endif

#define C(row,dim,j) ((row>>(j))&1)


typedef uint8_t  ind_t;
typedef uint64_t state_t;
typedef uint64_t vec_t;

static INLINE int scalar_product(vec_t a, vec_t b)
{
    return __builtin_parityl(a & b);
}

/* cache holds possible values for bit-sequences of orientable RBMs */
void populate_cache(vec_t **cache, size_t *size, const int dim);

/* allocate memory for RBM matrix */
vec_t *init(const ind_t dim);

/* --- get maximal possible state of the oriented RBM matrix --- */
static INLINE state_t get_max_state(ind_t dim)
{
    state_t max_state = 1;
    for (ind_t c = 0; c < dim - 1; ++c) {
        max_state <<= (dim - c - 2);
    }
    return max_state - 1;
}

/* set values of rows of RBM matrix dpending on state */
static INLINE void matrix_by_state(vec_t *mat, const vec_t *cache, const state_t state, const ind_t dim)
{
    state_t mask;
    for (int i=1, c=1, j=dim-3; j>=0; c+=i++, j--) {
        mask = ((state_t)1<<i) - 1;
        mat[j] = cache[(state>>(c-1))&mask];
    }
}

static inline bool is_upper_triangular(const vec_t *mat, const ind_t dim)
{
    vec_t row_mask = 1;
    for (ind_t i=1; i<dim; ++i) {
        row_mask |= (1u<<i);
        if ( mat[i] & row_mask ) {
            return false;
        }
    }
    return true;
}

/* main workers of this app */
bool is_spinc(const vec_t *mat, const ind_t dim);

static inline int row_sum(vec_t r)
{
#if VEC_T_SIZE == 64
    return __builtin_popcountl(r);
#else
    return __builtin_popcount(r);
#endif
}

static INLINE bool is_orientable(const vec_t *mat, const ind_t dim)
{
    for (ind_t j = 0; j < dim; ++j) {
        if (row_sum(mat[j]) & 1) {
            return false;
        }
    }
    return true;
}

bool is_spin(const vec_t *mat, const ind_t dim);

void print_mat(const vec_t *mat, const ind_t dim);

void swap_rows_and_cols(vec_t *src, vec_t *dst, ind_t dim, ind_t r1, ind_t r2);

void conditional_add_col(const vec_t *src, vec_t *dst, ind_t dim, ind_t k);

bool conditional_add_row(const vec_t *src, vec_t *dst, ind_t dim, ind_t l, ind_t m);

int matrix_weight(const vec_t *mat, const ind_t dim);

/* to delete - not needed anymore */

extern const ind_t mat_characters[];

void encode_matrix(const vec_t *mat, const ind_t dim, char *buffer);

void decode_matrix(vec_t *mat, const ind_t dim, const char *buffer);

char* matrix_to_digraph6(const vec_t *mat, int dim);

int digraph6_to_matrix(const char *digraph6_str, vec_t *matrix, int n);

#endif /* BOTT_H */
