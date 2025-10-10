#include <stdio.h>
#include "dag.h"

int main(void) {
    vec_t mat[64];
    int n;
    char s[2048];
    int all_even_out_degree;

    // Main loop to read one line (one graph) at a time
    while (fgets(s, sizeof(s), stdin) != NULL) {

        // 1. Determine the number of vertices (n) from the graph6 string.
        n = graphsize(s);
        if (n == 0) {
            continue; // Skip invalid or empty lines
        }
        matrix_from_d6(s, mat, n);
        // --- Filtering Logic ---
        all_even_out_degree = 1;
        for (int i = 0; i < n; ++i) {
            int out_degree = row_sum(mat[i]);
            if ((out_degree & 1) == 1) {
                all_even_out_degree = 0;
                break;
            }
        }

        if (all_even_out_degree) {
            printf("%s", s);
        }
    }
    return 0;
}
