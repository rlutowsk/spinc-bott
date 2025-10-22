#include "bott.h"
#include <nauty.h>
#include <assert.h>

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
    vec_t *mat = (vec_t*)malloc(sizeof(vec_t) * dim);
    memset(mat, 0, sizeof(vec_t) * dim);
    return mat;
}

void matrix_by_state(vec_t *mat, const vec_t *cache, const state_t state, const ind_t dim)
{
    state_t mask;
    for (int i=1, c=1, j=dim-3; j>=0; c+=i++, j--) {
        mask = ((state_t)1<<i) - 1;
        mat[j] = cache[(state>>(c-1))&mask];
    }
}

void print_mat(const vec_t *mat, const ind_t dim)
{
    static const char c[] = { '.', '1' };

    for (ind_t i=0; i<dim; i++) {
        for (ind_t j=0; j<dim; j++) {
            printf("%c", c[C(mat[i],dim,j)]);
        }
        printf("\n");
    }
}

static inline bool equal_cols(const vec_t *mat, const ind_t dim, const ind_t i, const ind_t j)
{
    vec_t mask = (1<<i) ^ (1<<j);
    for (ind_t r=0; r<dim; r++) {
        if (scalar_product(mat[r],mask)) {
            return false;
        }
    }
    return true;
}

bool is_upper_triangular(const vec_t *mat, const ind_t dim)
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
            if (e && aij!=C(mat[i],dim,j)) {
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
                if (scalar_product(mat[i], mat[j]) != C(mat[i],dim,j)) {
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

//************************************************************************************************
/* this is to delete, since we will use d6 format from dag.c file */

const vec_t row_masks[]      = { 0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095 };
const ind_t row_characters[] = { 0, 1, 1, 1,  1,  2,  2,   2,   2,   3,    3,    3,    3 };
const ind_t mat_characters[] = { 0, 1, 2, 3,  4, 10, 12,  14,  16,  27,   30,   33,   36 };

void encode_matrix(const vec_t *mat, const ind_t dim, char *buffer)
{
    for (ind_t i=0; i<dim; i++) {
        snprintf(buffer+i*row_characters[dim], row_characters[dim]+1, "%0*lx", row_characters[dim], mat[i]&row_masks[dim]);
    }
}

char* matrix_to_digraph6(const vec_t *mat, int dim) {
    // This function encodes a Bott matrix into a digraph6 string.

    // Ensure the output buffer is large enough
    // for the maximum size of dim (which is 12).
    // The maximum size is 1 (header) + 24 (body) + 1 (null terminator) = 26.
    // So a buffer of size 100 is more than sufficient.
    static char digraph_string[100];

    // We assume dim >= 0 and dim < 12, so no need for extensive checks.
    if (dim == 0) {
        digraph_string[0] = '@'; // ASCII 64, for n=0
        digraph_string[1] = '\0';
        return digraph_string;
    }

    // Number of bits in the adjacency matrix
    int num_bits = dim * dim;
    // Number of 6-bit chunks needed for the body
    int num_bytes_body = (num_bits + 5) / 6;
    int total_length = 1 + 1 + num_bytes_body + 1; // '&' + 1 for header, 1 for null terminator

    // Step 1: Encode the number of vertices (dim). Since dim <= 12,
    // for GAP digraphs start with '&' (ASCII 38).
    digraph_string[0] = '&'; // ASCII 38, for n <= 62
    // Next character encodes the number of vertices.
    digraph_string[1] = dim + 63;

    // Step 2: Encode the bit stream from the matrix directly into the string
    int bit_count = 0;
    int output_index = 2; // Start after '&' and header
    unsigned char current_byte = 0;

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            // Get the current bit
            int bit = C(mat[i], dim, j);

            // Pack the bit into the current 6-bit chunk
            current_byte |= (bit << (5 - (bit_count % 6)));
            bit_count++;

            // If a 6-bit chunk is complete, add it to the string
            if (bit_count % 6 == 0) {
                digraph_string[output_index] = current_byte + 63;
                output_index++;
                current_byte = 0;
            }
        }
    }

    // Pad the last chunk with zeros if it's not full and add to the string
    if (bit_count % 6 != 0) {
        digraph_string[output_index] = current_byte + 63;
        output_index++;
    }

    digraph_string[total_length - 1] = '\0';
    return digraph_string;
}

/**
 * @brief Decodes a single character from the digraph6 format.
 *
 * This helper function converts a character from the base-64 digraph6 alphabet
 * back to its 6-bit integer value.
 *
 * @param c The character to decode.
 * @return The 6-bit value.
 */
static inline unsigned char d6_char_to_int(char c) {
    if (c >= '?' && c <= '~') {
        return c - '?';
    } else if (c >= '!' && c <= '`') {
        return c - '!';
    }
    return 0; // Should not happen with valid digraph6 strings
}

/**
 * @brief Converts a digraph6 string into an adjacency matrix.
 *
 * This function manually parses a digraph6 string and populates a pre-allocated
 * matrix with the adjacency information. The matrix is represented as an array
 * of vec_t, where each vec_t corresponds to a row of the adjacency matrix.
 *
 * @param digraph6_str The input digraph6 string.
 * @param matrix A pre-allocated array of vec_t to store the adjacency matrix.
 * It must be of size n.
 * @param n The number of vertices in the graph.
 * @return 0 on success, -1 on failure (e.g., if the string is too short).
 */
int digraph6_to_matrix(const char *digraph6_str, vec_t *matrix, int n) {
    if (n >= 63) {
        // We can't use uint64_t for rows if n > 64, so we return an error.
        return -1;
    }

    // Initialize the matrix with zeros
    memset(matrix, 0, n * sizeof(vec_t));

    // Determine the start of the data stream, skipping the header.
    // The d6 format is variable-length for the number of vertices.
    const char *data_ptr = digraph6_str + 2;


    // This is the corrected logic. We use a bit index to track our position
    // in the bitstream and iterate through each character of the d6 string.
    int current_bit = 5;
    unsigned char current_char_value = d6_char_to_int(*data_ptr);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // if (i == j) {
            //     continue; // Skip the diagonal
            // }
            
            int bit = (current_char_value >> current_bit) & 1;
            
            if (bit) {
                matrix[i] |= (1ULL << j);
            }
            
            current_bit--;
            if (current_bit < 0) {
                data_ptr++;
                current_char_value = d6_char_to_int(*data_ptr);
                current_bit = 5;
            }
        }
    }

    return 0;
}

void decode_matrix(vec_t *mat, const ind_t dim, const char *buffer)
{
    if (dim == 0) {
        return;
    }

    const size_t row_chars = row_characters[dim];
    char temp_buffer[row_chars + 1];

    for (ind_t i = 0; i < dim; i++) {
        memcpy(temp_buffer, buffer + i * row_chars, row_chars);
        temp_buffer[row_chars] = '\0';

        mat[i] = strtol(temp_buffer, NULL, 16);
    }
}
