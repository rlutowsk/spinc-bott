#include "common.h"
#include "dag.h"
#include "bucket.h"
#include "adjpack11.h"
#include "parse_scaled.h"
#include "common.h"
#include "tlsbuf.h"

#include <omp.h>

void help(const char *progname) {
    if (progname == NULL) progname = "uniqueg";
    printf("Usage: %s [options]\n", progname);
    printf("Read d6-packed graphs from stdin and write unique codes to stdout.\n\n");
    printf("Options:\n");
    printf("  -c           Canonicalize graphs before checking uniqueness.\n");
    printf("  -v           Increase verbosity (can be repeated).\n");
    printf("  -j NUM       Number of OpenMP threads (default: %d).\n", omp_get_max_threads());
    printf("  -n NUM       Number of hash shards (default: 256).\n");
    printf("  -l SIZE      Batch size (number of lines) to process at once (default: 100000).\n");
    printf("               SIZE accepts suffixes like k, M, G or binary Ki, Mi (examples: 500k, 2M, 1G, 2Mi).\n");
    printf("  -i FILE      Read input from FILE instead of stdin.\n");
    printf("  -o FILE      Write output to FILE instead of stdout.\n");
    printf("  -h           Show this help message and exit.\n\n");
    printf("Examples:\n");
    printf("  cat graphs.d6 | %s -c -n 512 > uniques.d6\n", progname);
    printf("  %s -j 8 -l 500k -i input.d6 -o uniques.d6\n", progname);
}

int main(int argc, char *argv[]) {

    // always zero the timer at start of main
    tic();
    // default settings
    size_t lines_capacity = 100000; 
    int num_shards = 256;
    int num_threads = omp_get_max_threads(); // default max threads
    int v = 0;
    bool canon = false;

    FILE *in = stdin, *out = stdout;

    int opt = 1;
    while ((opt = getopt(argc, argv, "l:n:j:vi:o:c")) != -1) {
        switch (opt) {
        case 'v':
            ++v;
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
        case 'c':
            canon = true;
            break;
        case 'h':
            help(argv[0]);
            exit(EXIT_SUCCESS);
        default: /* '?' */
            help(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    // variables
    verbosity_level = v;

    size_t line_count = 0;
    char buffer[MAXLINE];

    printlog(1, "Using %d threads, %d shards, batch size %zu", num_threads, num_shards, lines_capacity);

    char **lines = malloc(lines_capacity * sizeof(char*));
    assert(lines != NULL);

#if !NAUTY_HAS_TLS
    fprintf(stderr, "Warning: nauty without TLS support, using single thread!\n");
    num_threads = 1;
#endif
    omp_set_num_threads(num_threads);

    if (!fgets(buffer, MAXLINE, in)) {
        fprintf(stderr, "error reading input\n");
        exit(1);
    }
    lines[line_count++] = strdup(buffer);

    int dim = graphsize(buffer);
    assert(dim > 0 && dim <= 11);

    // Global dedup across orbits by CANONICAL key
    GHashBucket *g_canonical_set =g_bucket_new_128(g_free, NULL, num_shards);

    size_t batch_num = 0;
    size_t num_of_reps = 0;

    while (true) {
        // read lines_capacity lines or until EOF
        while (line_count < lines_capacity && fgets(buffer, MAXLINE, in)) {
            lines[line_count++] = strdup(buffer);
        }
        if (line_count == 0) break;

        printlog(2, "Processing batch %zu; reps found: %zu.", ++batch_num, ATOMIC_GET(num_of_reps));

        #pragma omp parallel
        {
            OutputBuffer thread_buffer;
            buffer_init(&thread_buffer, 1000, out);

            init_nauty_data(dim);

            size_t local_reps = 0;

            char *line;
            key128_t key;
            unsigned n = 0;

            #pragma omp for schedule(dynamic, 100)
            for (size_t i = 0; i < line_count; ++i) {
                line = lines[i];

                d6pack_decode(line, &key, &n);
                if (canon){
                    vec_t m[11], c[11]; // max n = 11
                    adjpack_to_matrix(&key, m, n);
                    // Canonical key
                    matrix_to_matrix_canon(m, n, c);
                    adjpack_from_matrix(c, n, &key);
                }
                if ( g_bucket_insert_copy128(g_canonical_set, &key) ) {
                    ++local_reps;
                    buffer_add(&thread_buffer, line);
                }
            }
            
            #pragma omp atomic
                num_of_reps += local_reps;

            free_nauty_data();

            buffer_flush(&thread_buffer);
            buffer_destroy(&thread_buffer);
        }
        for (size_t i = 0; i < line_count; ++i) free(lines[i]);
        line_count = 0; // reset for next batch
    }

    free(lines);

    printlog(1, "Done. Found %lu representatives", num_of_reps); //g_bucket_size(g_canonical_set));

    g_bucket_destroy(g_canonical_set);

    if (in != stdin) fclose(in);
    if (out != stdout) fclose(out);

    return 0;
}
