#include <stdio.h>
#include "dag.h"

int main(void) {
    vec_t mat[64];
    int n;
    char s[2048];

    // Main loop to read one line (one graph) at a time
    while (fgets(s, sizeof(s), stdin) != NULL) {

        // 1. Determine the number of vertices (n) from the graph6 string.
        n = graphsize(s);
        if (n == 0) {
            continue; // Skip invalid or empty lines
        }
        matrix_from_d6(s, mat, n);

        if (is_upper_triangular(mat, n)) {
            printf("%s", s);
        }
    }
    return 0;
}
