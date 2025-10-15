#ifndef BOTT_H
#define BOTT_H

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <omp.h>

/* ensure we are in 64-bit environment */
#include <bits/wordsize.h>
#if __WORDSIZE != 64
    #error "compile only in 64-bit mode"
#endif

#ifndef SWAP
    #define SWAP(T, a, b) do { T SWAP_TMP = a; a = b; b = SWAP_TMP; } while (0)
#endif

#define C(row,dim,j) ((row>>(j))&1)

#include <stdint.h>

typedef uint8_t  ind_t;
typedef uint64_t state_t;

#if VEC_T_SIZE == 64
    typedef uint64_t vec_t;
#elif VEC_T_SIZE == 32
    typedef uint32_t vec_t;
#else
    #error "wrong ind_t type"
#endif

#if VEC_T_SIZE == 64
    #define scalar_product(a,b) (vec_t)__builtin_parityl((a)&(b))
    //#define row_sum __builtin_popcountl
#else
    #define scalar_product(a,b) (vec_t)__builtin_parity((a)&(b))
    //#define row_sum __builtin_popcount
#endif

/* cache holds possible values for bit-sequences of orientable RBMs */
void populate_cache(vec_t **cache, size_t *size, const ind_t dim);

/* allocate memory for RBM matrix */
vec_t *init(const ind_t dim);

/* set values of rows of RBM matrix dpending on state */
void matrix_by_state(vec_t *mat, const vec_t *cache, const state_t state, const ind_t dim);

int is_upper_triangular(const vec_t *mat, const ind_t dim);

/* main workers of this app */
size_t is_spinc(const vec_t *mat, const ind_t dim);

static inline int row_sum(vec_t r)
{
#if VEC_T_SIZE == 64
    return __builtin_popcountl(r);
#else
    return __builtin_popcount(r);
#endif
}

int is_spin(const vec_t *mat, const ind_t dim);

void print_mat(const vec_t *mat, const ind_t dim);

extern const ind_t mat_characters[];

void encode_matrix(const vec_t *mat, const ind_t dim, char *buffer);

void decode_matrix(vec_t *mat, const ind_t dim, const char *buffer);

char* matrix_to_digraph6(const vec_t *mat, int dim);

int digraph6_to_matrix(const char *digraph6_str, vec_t *matrix, int n);

void swap_rows_and_cols(vec_t *src, vec_t *dst, ind_t dim, ind_t r1, ind_t r2);

void conditional_add_col(const vec_t *src, vec_t *dst, ind_t dim, ind_t k);

int conditional_add_row(const vec_t *src, vec_t *dst, ind_t dim, ind_t l, ind_t m);

int matrix_weight(const vec_t *mat, const ind_t dim);

#endif /* BOTT_H */
