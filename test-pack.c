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
    vec_t mat[11], out[11]; // max n = 11
    unsigned n;
    size_t match = 0, unmatch = 0;

    FILE *in = fopen("test.d6", "r");
    if (in == NULL) {
        in = stdin;
        printf("=== [pack] testing packing implementations for digraph6 strings from stdin ===\n");
    } else {
        printf("=== [pack] testing packing implementations for digraph6 strings from file 'test.d6' ===\n");
    }
    fflush(stdout);

    while (fgets(line, sizeof(line), in)) {

        // Decode digraph6 string into key128_t
        if (!d6pack_decode(line, &key_from_d6, &n)) {
            line[strcspn(line, "\r\n")] = 0;
            fprintf(stderr, "    Invalid digraph6 string: %s\n", line);
            continue;
        }

        // Convert digraph6 string to adjacency matrix
        if (matrix_from_d6(line, mat, n)<0) {
            line[strcspn(line, "\r\n")] = 0;
            fprintf(stderr, "    matrix_from_d6 failed for: %s\n", line);
            continue;
        }

        // Pack matrix into key128_t
        adjpack_from_matrix(mat, n, &key_from_mat);

        // Compare both keys
        if (memcmp(&key_from_d6, &key_from_mat, sizeof(key128_t)) != 0) {
            printf("    Mismatch for input: %s\n", line);
            ++unmatch;
            print_key_bits(&key_from_d6, "    d6:  ");
            print_key_bits(&key_from_mat,"    mat: ");
            exit(1);
        } else {
            ++match;
        }
        adjpack_to_matrix(&key_from_mat, out, n);
        if (memcmp(mat, out, n*sizeof(vec_t)) != 0) {
            fprintf(stderr, "    packing - unpacking error\n");
            exit(1);
        }
    }

    printf("=== [pack] all tests passed. ===\n");

    if (in != stdin) {
        fclose(in);
    }
    return 0;
}
