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
	char d6[128];
    key128_t key_from_d6, key_from_mat;
    vec_t mat[11], out[11]; // max n = 11
    unsigned n;
    size_t match = 0, unmatch = 0;

    if (fgets(line, sizeof(line), stdin) == NULL) {
        fprintf(stderr, "error reading first line, quitting...\n");
        exit(1);
    }
    n = graphsize(line);
    init_nauty_data(n);

    while (fgets(line, sizeof(line), stdin)) {
        // Remove trailing newline
        line[strcspn(line, "\n")] = 0;

		d6_to_d6_canon(line, d6);
        // Decode digraph6 string into key128_t
        if (!d6pack_decode(d6, &key_from_d6, &n)) {
            fprintf(stderr, "Invalid digraph6 string: %s\n", line);
            continue;
        }

        // Convert digraph6 string to adjacency matrix
        if (matrix_from_d6(line, mat, n)<0) {
            fprintf(stderr, "matrix_from_d6 failed for: %s\n", line);
            continue;
        }
		matrix_to_matrix_canon(mat, n, out);
        // Pack matrix into key128_t
        adjpack_from_matrix(out, n, &key_from_mat);

        // Compare both keys
        if (memcmp(&key_from_d6, &key_from_mat, sizeof(key128_t)) != 0) {
            printf("Mismatch for input: %s\n", line);
            ++unmatch;
            print_key_bits(&key_from_d6, "d6:  ");
            print_key_bits(&key_from_mat,"mat: ");
            exit(1);
        } else {
            ++match;
        }
    }
    free_nauty_data();
    printf("Matched: %lu. Unmatched: %lu.\n", match, unmatch);

    return 0;
}
