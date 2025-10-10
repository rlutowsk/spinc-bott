#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "dag.h"
#include "adjpack11.h"
#include "d6pack11.h"  // for d6pack_decode

void print_key_bits(const key128_t *k, const char *label) {
    printf("%s: ", label);
    for (int i = 0; i < 16; ++i) {
        for (int b = 7; b >= 0; --b) {
            putchar((k->b[i] >> b) & 1 ? '1' : '0');
        }
    }
    putchar('\n');
}

int main(void) {
    char line[128];
    key128_t key_from_d6, key_from_mat;
    vec_t mat[11]; // max n = 11
    unsigned n;
    size_t match = 0, unmatch = 0;

    while (fgets(line, sizeof(line), stdin)) {
        // Remove trailing newline
        line[strcspn(line, "\n")] = 0;

        // Decode digraph6 string into key128_t
        if (!d6pack_decode(line, &key_from_d6, &n)) {
            fprintf(stderr, "Invalid digraph6 string: %s\n", line);
            continue;
        }

        // Convert digraph6 string to adjacency matrix
        if (matrix_from_d6(line, mat, n)<0) {
            fprintf(stderr, "matrix_from_d6 failed for: %s\n", line);
            continue;
        }

        // Pack matrix into key128_t
        adjpack_from_matrix(mat, n, &key_from_mat);

        // Compare both keys
        if (memcmp(&key_from_d6, &key_from_mat, sizeof(key128_t)) != 0) {
            printf("Mismatch for input: %s\n", line);
            ++unmatch;
        } else {
            ++match;
        }
    }

    printf("Matched: %lu. Unmatched: %lu.\n", match, unmatch);

    return 0;
}