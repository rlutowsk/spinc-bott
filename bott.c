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

void print_mat(const vec_t *mat, const ind_t dim)
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

const vec_t row_masks[]      = { 0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095 };
const ind_t row_characters[] = { 0, 1, 1, 1,  1,  2,  2,   2,   2,   3,    3,    3,    3 };
const ind_t mat_characters[] = { 0, 1, 2, 3,  4, 10, 12,  14,  16,  27,   30,   33,   36 };

void encode_matrix(const vec_t *mat, const ind_t dim, char *buffer)
{
    ind_t i;
    for (i=0; i<dim; i++) {
        snprintf(buffer+i*row_characters[dim], row_characters[dim]+1, "%0*lx", row_characters[dim], mat[i]&row_masks[dim]);
    }
}

void decode_matrix(vec_t *mat, const ind_t dim, const char *buffer)
{
    ind_t i;
    for (i = 0; i < dim; i++) {
        // Obliczanie przesunięcia w buforze
        size_t offset = i * row_characters[dim];

        // Tymczasowy bufor do podciągu
        char temp_buffer[row_characters[dim] + 1];
        strncpy(temp_buffer, buffer + offset, row_characters[dim]);
        temp_buffer[row_characters[dim]] = '\0';

        // Konwersja z szesnastkowego ciągu na liczbę
        mat[i] = strtol(temp_buffer, NULL, 16);
    }
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
void conditional_add_col(vec_t *src, vec_t *dst, ind_t dim, ind_t k)
{
    vec_t row_k = src[k];

    for (ind_t i=0; i<dim; i++) {
        vec_t cond = -(vec_t)(((src[i] >> k) & 1));
        dst[i] = src[i] ^ (row_k & cond);
    }
}

/* Op3 */
void conditional_add_row(vec_t *src, vec_t *dst, ind_t dim, ind_t l, ind_t m)
{
    vec_t row_l = src[l];

    for (ind_t i=0; i<dim; i++) {
        dst[i] = src[i];
    }
    vec_t cond = -(vec_t)((l!=m) && equal_cols(src, dim, l, m));
    dst[m] = src[m] ^ (row_l & cond);
}