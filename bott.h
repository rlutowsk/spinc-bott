#ifndef BOTT_H
#define BOTT_H

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

#ifndef SWAP
# define SWAP(T, a, b) do { T SWAP_TMP = a; a = b; b = SWAP_TMP; } while (0)
#endif

#define C(row,dim,j) ((row>>(dim-j-1))&1)

#define VEC_T_SIZE 64

#include <stdint.h>

typedef uint8_t  ind_t;
typedef uint64_t state_t;

#if VEC_T_SIZE == 64
    typedef uint64_t vec_t;
#elif VEC_T_SIZE == 32
    typedef uint32_t vec_t;
#else
#   error "wrong ind_t type"
#endif

#if VEC_T_SIZE == 64
#define scalar_product(a,b) __builtin_parityl((a)&(b))
//#define row_sum __builtin_popcountl
#else
#define scalar_product(a,b) __builtin_parity((a)&(b))
//#define row_sum __builtin_popcount
#endif

/* cache holds possible values for bit-sequences of orientable RBMs */
void populate_cache(vec_t **cache, size_t *size, const ind_t dim);

/* allocate memory for RBM matrix */
vec_t *init(const ind_t dim);

/* set values of rows of RBM matrix dpending on state */
void set(vec_t *mat, const vec_t *cache, const state_t state, const ind_t dim);

/* main workers of this app */
size_t is_spinc(const vec_t *mat, const ind_t dim);
size_t is_spin(const vec_t *mat, const ind_t dim);

#endif
