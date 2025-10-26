#include "common.h"
#include "dag.h"
#include "adjpack11.h"

int main(void) {
    vec_t mat[11];
    unsigned n = 0;
    char s[MAXLINE];
    key128_t k;

    // Main loop to read one line (one graph) at a time
    while (fgets(s, sizeof(s), stdin) != NULL) {
        remove_newline(s);
        // 1. Determine the number of vertices (n) from the graph6 string.
        if (!d6pack_decode(s, &k, &n)) {
            continue; // Skip invalid or empty lines
        }
        adjpack_to_matrix(&k, mat, n);

        if (!is_upper_triangular(mat,n)) {
            fprintf(stderr, "all input matrices must be upper triangular, corrupted input, quitting...\n");
            exit(1);
        }

        if (is_spinc(mat, n)) {
            printf("%s\n", s);
        }
    }
    return 0;
}
