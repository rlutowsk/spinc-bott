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


int matrix_from_d6_1(char *s, vec_t *mat, ind_t dim)
{
    char *p;
    int n,i,j,k,x;

    if (s[0] != '&') {
        return -1;
    }

    n = graphsize(s);
    
    if (n <= 0 || n > dim) {
        return -1;
    }

    memset(mat, 0, dim * sizeof(vec_t)); // Initialize the matrix with zeros

    x = 0;
    p = s + 1 + SIZELEN(n);
    // p = s + SIZELEN(n);

    k = 1;
    for (i = 0; i < n; ++i)
    {
        // mat[i] = 0;
        for (j = 0; j < n; ++j)
        {
            if (--k == 0)
            {
                k = 6;
                x = *(p++) - BIAS6;
            }
            
            // if ((x & TOPBIT6))
            // {
            //     mat[i] |= (1ULL << j);
            // }
            
            mat[i] |= (vec_t)((x & TOPBIT6)!=0) << j;
            x <<= 1;
        }
    }
    return 0;
}

int matrix_from_d6_2(char *s, vec_t * __restrict mat, ind_t dim)
{
    key128_t k;
    unsigned n = 0;
    d6pack_decode(s, &k, &n);
    if (n == 0 || (ind_t)n > dim) return -1;
    adjpack_to_matrix(&k, mat, n);
    return 0;
}

int main(void) {
    char line[MAXLINE];
    vec_t mat_1[11], mat_2[11]; // max n = 11
    size_t match = 0, unmatch = 0;
    uint64_t t1, t2,
        t_matrix_from_d6_1 = 0,
        t_matrix_from_d6_2 = 0;
    int rval;

    unsigned n = 1;
    init_nauty_data(n);

    FILE *in = fopen("test.d6", "r");
    if (in == NULL) {
        in = stdin;
        printf("=== [matrix_from_d6] testing implementations for digraph6 strings from stdin ===\n");
    } else {
        printf("=== [matrix_from_d6] testing implementations for digraph6 strings from file 'test.d6' ===\n");
    }
    fflush(stdout);
    while (fgets(line, sizeof(line), in)) {
        if (n != graphsize(line)) {
            n = graphsize(line);
            free_nauty_data();
            init_nauty_data(n);
        }

        // Convert digraph6 string to adjacency matrix
        t1 = ns_now_monotonic();
        rval = matrix_from_d6_1(line, mat_1, n);
        t2 = ns_now_monotonic();
        t_matrix_from_d6_1 += (t2 - t1);
        if (rval < 0) {
            fprintf(stderr, "    matrix_from_d6_1 failed for: %s\n", line);
            exit(1);
        }

        t1 = ns_now_monotonic();
        rval = matrix_from_d6_2(line, mat_2, n);
        t2 = ns_now_monotonic();
        t_matrix_from_d6_2 += (t2 - t1);
        if (rval < 0) {
            fprintf(stderr, "    matrix_from_d6_2 failed for: %s\n", line);
            exit(1);
        }
        // Compare both matrices
        if (memcmp(mat_1, mat_2, n * sizeof(vec_t)) != 0) {
            printf("    Mismatch for input: %s\n", line);
            ++unmatch;
            exit(1);
        } else {
            ++match;
        }
    }
    free_nauty_data();
    if (in != stdin) {
        fclose(in);
    }
    printf("=== [matrix_from_d6] timings (ms):\n");
    printf("    matrix_from_d6_1:  %.3f\n", t_matrix_from_d6_1/1e6);
    printf("    matrix_from_d6_2:  %.3f\n", t_matrix_from_d6_2/1e6);
    printf("=== [matrix_from_d6] all tests completed ===\n");
   
    return 0;
}
