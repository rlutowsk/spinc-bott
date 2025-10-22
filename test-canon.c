#define _GNU_SOURCE
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "common.h"
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
    size_t match = 0, unmatch = 0;
    uint64_t t1, t2,
        t_d6_to_canon = 0,
        t_matrix_from_d6 = 0,
        t_matrix_to_matrix_canon = 0,
        t_adjpack_from_matrix = 0,
        t_d6pack_decode = 0;
    int rval;


    unsigned n = 1;
    init_nauty_data(n);
    atexit(free_nauty_data);

    FILE *in = fopen("test.d6", "r");
    if (in == NULL) {
        in = stdin;
        printf("=== [canon] testing digraph6 strings from stdin ===\n");
    } else {
        printf("=== [canon] testing digraph6 strings from file 'test.d6' ===\n");
    }
    fflush(stdout);

    while (fgets(line, sizeof(line), in)) {
        if (n != graphsize(line)) {
            n = graphsize(line);
            free_nauty_data();
            init_nauty_data(n);
        }

        t1 = ns_now_monotonic();
		d6_to_d6_canon(line, d6);
        t2 = ns_now_monotonic();
        t_d6_to_canon += (t2 - t1);

        // Decode digraph6 string into key128_t
        t1 = ns_now_monotonic();
        rval = d6pack_decode(d6, &key_from_d6, &n);
        t2 = ns_now_monotonic();
        t_d6pack_decode += (t2 - t1);

        if (!rval) {
            fprintf(stderr, "    Invalid digraph6 string: %s\n", line);
            exit(1);
        }

        // Convert digraph6 string to adjacency matrix
        t1 = ns_now_monotonic();
        rval = matrix_from_d6(line, mat, n);
        t2 = ns_now_monotonic();
        t_matrix_from_d6 += (t2 - t1);
        if (rval < 0) {
            fprintf(stderr, "    matrix_from_d6 failed for: %s\n", line);
            exit(1);
        }

        t1 = ns_now_monotonic();
		matrix_to_matrix_canon(mat, n, out);
        t2 = ns_now_monotonic();
        t_matrix_to_matrix_canon += (t2 - t1);
        // Pack matrix into key128_t
        t1 = ns_now_monotonic();
        adjpack_from_matrix(out, n, &key_from_mat);
        t2 = ns_now_monotonic();
        t_adjpack_from_matrix += (t2 - t1);
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
    }
    free_nauty_data();

    if (in != stdin) {
        fclose(in);
    }

    printf("=== [canon] timings (ns) ===\n");
    printf("    d6_to_d6_canon:        %15lu\n", t_d6_to_canon);
    printf("    d6pack_decode:         %15lu\n", t_d6pack_decode);
    printf("    matrix_from_d6:        %15lu\n", t_matrix_from_d6);
    printf("    matrix_to_matrix_canon:%15lu\n", t_matrix_to_matrix_canon);
    printf("    adjpack_from_matrix:   %15lu\n", t_adjpack_from_matrix);
    printf("=== [canon] all tests passed ===\n");
    return 0;
}
