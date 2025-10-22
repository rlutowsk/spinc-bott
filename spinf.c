#include <stdio.h>
#include "dag.h"
#include "bott.h"
#include "adjpack11.h"

int main(void) {
    vec_t mat[11];
    unsigned n;
    char s[128];
    key128_t k;

    // Main loop to read one line (one graph) at a time
    while (fgets(s, sizeof(s), stdin) != NULL) {
        s[strcspn(s, "\r\n")] = 0;

        // 1. Write d6 code to key128_t k and determine the number of vertices (n).
        if (!d6pack_decode(s, &k, &n)) {
            continue; // Skip invalid or empty lines
        }

        // 2. Convert the key128_t k to an adjacency matrix in vec_t mat.
        adjpack_to_matrix(&k, mat, n);

        // 3. Check if the matrix is upper triangular.
        if (!is_upper_triangular(mat,n)) {
            fprintf(stderr, "all input matrices must be upper triangular, corrupted input, quitting...\n");
            exit(1);
        }

        // 4. Check if the matrix represents a spin real Bott manifold.
        if (is_spin(mat, n)) {
            printf("%s\n", s);
        }
    }
    return 0;
}
