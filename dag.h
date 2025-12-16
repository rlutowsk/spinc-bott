#pragma once

#include <nauty.h>
#include <gtools.h>
#include <stdbool.h>

#include "d6pack11.h"

/**
 * @brief Converts a nauty graph to a digraph6 (d6) string representation.
 *
 * This function encodes the given nauty graph into the digraph6 format,
 * including the null terminator. The resulting string is stored in a
 * dynamically allocated buffer, which should not be freed by the caller
 * (the buffer is managed internally, likely as a global variable).
 *
 * @param g Pointer to the nauty graph structure.
 * @param m The number of setwords per row (as required by nauty).
 * @param n The number of vertices in the graph.
 * @return Pointer to a null-terminated string containing the digraph6 encoding of the graph.
 */
 char* graph_to_d6(graph *g, int m, int n, char *dag_gcode);

void init_nauty_data(int n);
void free_nauty_data(void);

/**
 * @brief Converts a binary matrix to a digraph6 (d6) string representation.
 *
 * The function encodes the given n x n adjacency matrix into a string using a
 * 6-bit encoding scheme, suitable for compact storage or transmission. The
 * resulting string begins with '&', followed by the encoded graph size, and
 * then the adjacency information packed into 6-bit characters.
 *
 * @param mat Pointer to the adjacency matrix (of type vec_t).
 * @param n   The number of vertices in the graph (matrix is n x n).
 * @return    Pointer to a statically allocated string containing the d6 encoding.
 *
 * Note: The returned pointer refers to a buffer managed by the function and may
 * be overwritten by subsequent calls. The caller should copy the result if it
 * needs to be preserved.
 */
char* matrix_to_d6(const vec_t *mat, int dim, char *dag_gcode);

/**
 * @brief Fills a matrix from a d6-encoded graph string.
 *
 * This function decodes a graph represented in the "d6" format (a compact
 * graph encoding) and fills the provided matrix with the adjacency information.
 * The matrix is assumed to be an array of `vec_t` values, where each entry
 * represents a row of the adjacency matrix as a bit vector.
 *
 * @param s    Pointer to the d6-encoded graph string. The string must start with '&'.
 * @param mat  Pointer to the matrix (array of `vec_t`) to be filled. The matrix
 *             should have space for at least `dim` rows.
 * @param dim  The dimension (number of rows/columns) of the matrix.
 *
 * @return 0 on success, -1 on error (e.g., invalid input string or dimension).
 *
 * The function expects the following:
 *   - The input string `s` is in the correct d6 format.
 *   - The matrix `mat` is pre-allocated and large enough to hold the result.
 *   - The function initializes the matrix to zero before filling it.
 */
int matrix_from_d6(char *s, vec_t *mat, ind_t dim);
/**
 * @brief Converts an adjacency matrix to its canonical digraph6 (d6) representation.
 *
 * Converts a given adjacency matrix representation of a directed graph (DAG) into
 * D6 string representation of its canonical form.
 *
 * @param mat Pointer to the adjacency matrix, where each row is represented as a vec_t.
 * @param n   The number of vertices in the graph (size of the matrix).
 * @return    A pointer to a string containing the canonical D6 representation of the graph.
 *
 * The function allocates and initializes a nauty graph structure from the input matrix,
 * skipping rows with no outgoing edges. It then generates the canonical form of the digraph
 * and returns its D6 string representation. The caller is responsible for freeing the returned string.
 */
char* matrix_to_d6_canon(const vec_t *mat, int dim, char *dag_gcode);

vec_t* matrix_to_matrix_canon(const vec_t *mat, int n, vec_t *out);

key128_t *matrix_to_key128_canon(const vec_t *mat, int n, key128_t *key);

char* d6_to_d6_canon(const char *src, char *dst);

key128_t* d6_to_key128_canon(const char *src, key128_t *key);

char *d6_to_d6_upper(char *src, char *dst);

int matrix_from_graph(graph *g, vec_t *mat);
int matrix_to_graph(graph *g, const vec_t *mat);

graph *get_dag_g(void);
