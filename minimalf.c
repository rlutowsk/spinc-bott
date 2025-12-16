#include <omp.h>

#include "bott.h"
#include "bucket.h"
#include "dag.h"
#include "adjpack11.h"
#include "parse_scaled.h"
#include "tlsbuf.h"

static void help(const char *progname) {
    fprintf(stderr,
        "Usage: %s [-j num_threads] [-n num_shards] [-l lines_capacity] [-i input] [-o output] [-v] [-u]\n"
        "  -j NUM   Number of threads (default: max available)\n"
        "  -n NUM   Number of hash shards (default: 256)\n"
        "  -l NUM   Input lines per batch (default: 100000, accepts k/M/G suffixes)\n"
        "  -i FILE  Input file (default: stdin)\n"
        "  -o FILE  Output file (default: stdout), use '-' for stdout also\n"
        "  -u       Output unique representatives only\n"
        "  -v       Increase verbosity (can be repeated)\n"
        "  -h       Show this help message\n",
        progname);
}

#define INITIAL_CAPACITY 1024

typedef struct {
    vec_t *rows;     // flat array of rows
    size_t len;      // number of matrices stored
    size_t capacity; // total number of matrices that can be stored
    size_t dim;      // rows per matrix
} MatArray;

static MatArray *matarray_create(size_t dim) {
    MatArray *arr = malloc(sizeof(MatArray));
    arr->dim = dim; arr->len = 0; arr->capacity = INITIAL_CAPACITY;
    arr->rows = malloc(sizeof(vec_t) * dim * arr->capacity);
    return arr;
}
static void matarray_free(MatArray *arr) {
    if (arr) { free(arr->rows); free(arr); }
}
static void matarray_append(MatArray *arr, const vec_t *mat) {
    if (arr->len >= arr->capacity) {
        arr->capacity *= 2;
        arr->rows = realloc(arr->rows, sizeof(vec_t) * arr->dim * arr->capacity);
        assert(arr->rows != NULL);
    }
    memcpy(&arr->rows[arr->len * arr->dim], mat, sizeof(vec_t) * arr->dim);
    arr->len++;
}
static vec_t *matarray_get(MatArray *arr, size_t index) {
    assert(index < arr->len);
    return &arr->rows[index * arr->dim];
}

/**
 * BFS over the forward closure under conditional_add_* (mod iso).
 * Returns false ASAP if any canonical neighbor key < canonical seed key.
 *
 * IMPORTANT: To match orbitg.c semantics, we:
 *  - enqueue NON-CANONICAL matrices produced by the operations,
 *  - canonicalize neighbors only to derive the key for visited/min check.
 */
static ind_t dim = 0;

static bool is_orbit_minimum(const key128_t *seedk) {
    // Per-thread "visited" set holds canonical keys
    FlatSet visited_set;
    flat_init(&visited_set, 1024);

    // Queue holds NON-CANONICAL matrices (like orbitg.c)
    MatArray *q = matarray_create(dim);

    vec_t initial_mat[dim],
          seed_can_mat[dim],
          aux[dim],
          canon_neighbor[dim];
    key128_t k;

    // Decode seed key -> matrix using adjpack (must match how seedk was built)
    adjpack_to_matrix(seedk, initial_mat, dim);

    // Canonical seed key for comparisons and visited
    matrix_to_matrix_canon(initial_mat, dim, seed_can_mat);

    key128_t seed_can_key;
    adjpack_from_matrix(seed_can_mat, dim, &seed_can_key);
    flat_insert(&visited_set, &seed_can_key);

    // Start BFS from the NON-CANONICAL seed (to match orbitg)
    matarray_append(q, initial_mat);

    bool is_min = true;

    for (size_t h = 0; h < q->len; ++h) {
        vec_t *cur;

        // Op2: column additions
        for (ind_t i = 0; i < dim; ++i) {
            cur = matarray_get(q, h);                   // re-fetch every time (safe if realloc)
            conditional_add_col(cur, aux, dim, i);      // aux: NON-CANON neighbor

            // Canonicalize aux to get the visited key
            // vec_t canon_neighbor[dim];
            matrix_to_matrix_canon(aux, dim, canon_neighbor);

            // key128_t k;
            adjpack_from_matrix(canon_neighbor, dim, &k);

            // Minimality test: neighbor < canonical seed?
            if (key128_lt(&k, &seed_can_key)) { is_min = false; goto cleanup; }

            // First time we see this canonical element? Enqueue the NON-CANON aux
            if (!flat_lookup(&visited_set, &k) && flat_insert(&visited_set, &k)) {
                matarray_append(q, aux);
            }
        }

        // Op3: row additions (both directions)
        for (ind_t i = 0; i < dim; ++i) {
            for (ind_t j = i + 1; j < dim; ++j) {
                cur = matarray_get(q, h);
                if (conditional_add_row(cur, aux, dim, i, j)) {
                    //vec_t canon_neighbor[dim];
                    matrix_to_matrix_canon(aux, dim, canon_neighbor);
                    //key128_t k;
                    adjpack_from_matrix(canon_neighbor, dim, &k);
                    if (key128_lt(&k, &seed_can_key)) { is_min = false; goto cleanup; }
                    if (!flat_lookup(&visited_set, &k) && flat_insert(&visited_set, &k)) {
                        matarray_append(q, aux);
                    }
                }

                cur = matarray_get(q, h);
                if (conditional_add_row(cur, aux, dim, j, i)) {
                    // vec_t canon_neighbor[dim];
                    matrix_to_matrix_canon(aux, dim, canon_neighbor);
                    // key128_t k;
                    adjpack_from_matrix(canon_neighbor, dim, &k);
                    if (key128_lt(&k, &seed_can_key)) { is_min = false; goto cleanup; }
                    if (!flat_lookup(&visited_set, &k) && flat_insert(&visited_set, &k)) {
                        matarray_append(q, aux);
                    }
                }
            }
        }
    }

cleanup:
    matarray_free(q);
    flat_free(&visited_set);
    return is_min;
}

int main(int argc, char *argv[]) {

    // always zero the timer at start of main
    tic();
    // default settings
    size_t lines_capacity = 100000;
    int num_shards = 256;
    int num_threads = omp_get_max_threads(); // default to max threads
    bool unique = false;

    double time_start = omp_get_wtime();

    FILE *in = stdin, *out = stdout;

    int opt = 1;
    while ((opt = getopt(argc, argv, "l:n:j:vi:o:u")) != -1) {
        switch (opt) {
        case 'v':
            increase_verbosity();
            break;
        case 'j':
            num_threads = (ind_t)atoi(optarg);
            break;
        case 'n':
            num_shards = (int)atoi(optarg);
            break;
        case 'l':
            size_t tmp = 0;
            if (!parse_scaled_size(optarg, &tmp) || tmp == 0) {
                fprintf(stderr, "Invalid -l value: %s (examples: 500k, 2M, 1G, 2Mi)\n", optarg);
                exit(EXIT_FAILURE);
            }
            lines_capacity = tmp;
            break;
        case 'i':
            in = fopen(optarg, "r");
            if (in == NULL) {
                perror("Error opening input file");
                exit(EXIT_FAILURE);
            }
            break;
        case 'o':
            out = fopen(optarg, "w");
            if (out == NULL) {
                perror("Error opening output file");
                exit(EXIT_FAILURE);
            }
            break;
        case 'u':
            unique = true;
            break;
        case 'h':
            help(argv[0]);
            exit(EXIT_SUCCESS);
        default: /* '?' */
            help(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    char *buffer = (char*)malloc(MAXLINE * lines_capacity * sizeof(char));
    assert(buffer != NULL);

    printlog(1, "Using %d threads, %d shards, batch size %zu", num_threads, num_shards, lines_capacity);

    char **lines = malloc(lines_capacity * sizeof(char*));
    assert(lines != NULL);
    for (size_t i = 0; i < lines_capacity; ++i) {
        lines[i] = buffer + i * MAXLINE;
        *lines[i] = '\0';
    }

#if !NAUTY_HAS_TLS
    fprintf(stderr, "Warning: nauty without TLS support, using single thread!\n");
    num_threads = 1;
#endif
    omp_set_num_threads(num_threads);

    if (!fgets(buffer, MAXLINE, in)) {
        fprintf(stderr, "error reading input\n");
        exit(1);
    }
    // lines[line_count++] = buffer;
    // lines[0] points to buffer anyway
    size_t line_count = 1;

    dim = graphsize(lines[0]);
    assert(dim > 0 && dim <= 11);

    // Global dedup across orbits by CANONICAL key
    GHashBucket *g_canonical_set = NULL;
    if (unique) {
        g_canonical_set = g_bucket_new_128(free, num_shards);
    }

    size_t batch_num = 0;
    size_t num_of_reps = 0;
    double read_time = 0, comp_time = 0;

    #pragma omp parallel
    {
        // Per thread buffer for output lines
        OutputBuffer thread_buffer;
        buffer_init(&thread_buffer, 5000, out);

        init_nauty_data(dim);

        for (;;) {
            double start = 0.0;

            // Single thread reads a batch
            #pragma omp single
            {
                start = omp_get_wtime();
                while (line_count < lines_capacity && fgets(lines[line_count], MAXLINE, in)) {
                    ++line_count;
                }
                read_time += omp_get_wtime() - start;

                if (line_count != 0) {
                    printlog(2, "Processing batch %zu; reps found: %zu.", ++batch_num, num_of_reps);
                }
            }

            // Barrier synchronizing — all see line_count
            #pragma omp barrier

            if (line_count == 0) {
                // End of data – exit per-thread loop
                break;
            }

            // Processing batch – each thread operates on its own buffer
            double comp_start = omp_get_wtime();

            size_t local_reps = 0;

            #pragma omp for schedule(dynamic,1000)
            for (size_t i = 0; i < line_count; ++i) {
                char *line = lines[i];

                key128_t seedk;
                d6_to_key128(line, &seedk);

                // Minimality test via forward orbit
                if (is_orbit_minimum(&seedk)) {
                    if (unique) {
                        vec_t m[dim], c[dim];
                        adjpack_to_matrix(&seedk, m, dim);
                        // Canonical key for global dedup
                        matrix_to_matrix_canon(m, dim, c);

                        key128_t seed_can_key;
                        adjpack_from_matrix(c, dim, &seed_can_key);

                        if (g_bucket_insert_copy128(g_canonical_set, &seed_can_key)) {
                            ++local_reps;
                            buffer_add(&thread_buffer, line);   // *** per-wątek buf
                        }
                    } else {
                        ++local_reps;
                        buffer_add(&thread_buffer, line);
                    }
                }
            }

            // Update global count once per batch (after reduction)
            #pragma omp atomic
            num_of_reps += local_reps;

            // *** All threads done computing
            #pragma omp barrier

            #pragma omp single
            {
                comp_time += omp_get_wtime() - comp_start;
                line_count = 0;
            }

            // Closure of batch (not always necessary, but keeps things tidy)
            #pragma omp barrier
        }

        free_nauty_data();
        buffer_flush(&thread_buffer);
        buffer_destroy(&thread_buffer);
    }

    free(lines);
    free(buffer);

    double time_end = omp_get_wtime();
    printlog(1, "Times. Reading: %.3fs. Computations: %.3fs. Ratio: %.4f. Total: %.3fs. Ratio: %.4f", read_time, comp_time, comp_time/(read_time+comp_time), time_end - time_start, comp_time/(time_end - time_start));
    printlog(1, "Done. Found %lu representatives", num_of_reps);

    if (g_canonical_set) {
        g_bucket_destroy(g_canonical_set);
    }

    if (in != stdin) fclose(in);
    if (out != stdout) fclose(out);

    return 0;
}
