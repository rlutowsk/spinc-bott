#include "bott.h"

mat_t col_masks[11];
mat_t row_masks[11];

void populate_cache(vec_t **cache, size_t *size, const ind_t dim)
{
    vec_t i, j, o, m;

    *size = 1<<(dim-2);
    if (*size<4) { *size = 4; }
    *cache = (vec_t*)malloc(VEC_T_SIZE*(*size));
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
    mat = (vec_t*)malloc(VEC_T_SIZE*dim);
    memset(mat, 0, VEC_T_SIZE*dim);
    return mat;
}

void new_init(ind_t dim)
{
    if (dim > 11) {
        dim = 11;
    }
    for (int i = 0; i < dim; i++) {
        mat_t mask = 0;
        for (int k = 0; k < dim; k++) {
            mask |= (mat_t)1 << (k * dim + dim - i - 1);
        }
        col_masks[i] = mask;
        row_masks[i] = ( ((mat_t)1 << dim) - 1 ) << i*dim;
    }
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

void new_set(mat_t *mat, const vec_t *cache, const state_t state, const ind_t dim)
{
    mat_t result = 0;
    ind_t num_bits = 0;
    ind_t row_idx = dim - 3;

    // Przechodzi przez kolejne bloki bitów w state
    for (ind_t i = 1; ; i++) {
        // Oblicza pozycję, z której należy odczytać dane
        ind_t shift = num_bits;
        
        // Zwiększa liczbę bitów dla kolejnego kroku
        num_bits += i;

        // Jeśli przesunięcie jest większe niż rozmiar state_t, kończy pętlę
        if (shift >= sizeof(state_t) * 8) {
            break;
        }

        // Tworzy maskę dla i bitów
        state_t mask = ((state_t)1 << i) - 1;
        
        // Wyodrębnia wartość, używając maski
        vec_t index_val = (state >> shift) & mask;
        
        // Odczytuje wiersz z cache
        mat_t row = cache[index_val];
        
        // Zapisuje wiersz na odpowiedniej pozycji w macierzy
        result |= row << (row_idx * dim);

        if (row_idx==0) {
            break;
        }
        // Zmniejsza indeks wiersza
        row_idx--;
    }
    
    *mat = result;
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
    vec_t mask = (1<<(dim-i-1)) ^ (1<<(dim-j-1));
    ind_t r;
    for (r=0; r<i; r++) {
        if (scalar_product(mat[r],mask)) {
            return 0;
        }
    }
    /* r = i; */
    for (; r<j; r++) {
        if (mat[r]&mask) {
            return 0;
        }
    }
    return 1;
}

static inline elm_t new_equal_cols(mat_t mat, const ind_t dim, const ind_t i, const ind_t j)
{
    return COL(mat, dim, i) == COL(mat, dim, j);
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

size_t new_is_spinc(const mat_t mat, const ind_t dim)
{
    ind_t i, j, e, z;
    elm_t aij;
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
            aij = new_equal_cols(mat, dim, i, j) ? 0 : scalar_product(ROW(mat, dim, i), ROW(mat, dim, j));
            if (z && aij) {
                z = 0;
            }
            if (e && aij!=ELM(mat, dim, i, j)) {
                e = 0;
            }
            if (e==0 && z==0) {
                return 0;
            }
        }
    }
    return 1;
}

size_t debug_is_spinc(const vec_t *old, const mat_t mat, const ind_t dim)
{
    ind_t i, j, e, z;
    elm_t aij;
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
            aij = new_equal_cols(mat, dim, i, j) ? 0 : scalar_product(ROW(mat, dim, i), ROW(mat, dim, j));
            if (z && aij) {
                z = 0;
            }
            if (e && aij!=ELM(mat, dim, i, j)) {
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

