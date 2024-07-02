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

typedef unsigned char ind_t;
typedef unsigned long vec_t;

/* cache holds possible values for bit-sequences of orientable RBMs */
void populate_cache(vec_t **cache, size_t *size, const ind_t dim);

/* allocate memory for RBM matrix */
vec_t *init(const ind_t dim);

/* set values of rows of RBM matrix dpending on state */
void set(vec_t *mat, const vec_t *cache, const vec_t state, const ind_t dim);

/* main workers of this app */
size_t is_spinc(const vec_t *mat, const ind_t dim);
size_t is_spin(const vec_t *mat, const ind_t dim);

#endif