#include "bott.h"
#include "common.h"
#include "dag.h"
#include "bucket.h"

#define MAXLINE 128

static ind_t dim = 0;

void help(const char *name)
{
    fprintf(stderr, "Usage: %s [-i input] [-o output] [-v] [-h] [-n shards] [-t interval] [-m cap]\n", name);
    fprintf(stderr, "  -i input    Input file (default: stdin)\n");
    fprintf(stderr, "  -o output   Output file (default: stdout)\n");
    fprintf(stderr, "  -v          Increase verbosity (can be used multiple times)\n");
    fprintf(stderr, "  -h          Show this help message and exit\n");
    fprintf(stderr, "  -n shards   Number of hash table shards (default: 1023)\n");
    fprintf(stderr, "  -t interval Print stats every <interval> seconds (default: 1)\n");
    fprintf(stderr, "  -m cap      Set peak memory cap in MB (calculates shards)\n");
}
/**
 * @brief Checks if the output of local calculation is in the code set.
 *
 *
 * @param b     Global bucket that stores the d6 codes.
 * @param q     Queue for storing orbit of the current code.
 * @param aux   Auxiliary matrix stores temporary output of the calculations.
 * @return void
 */
static inline void add_code(GHashBucket *b, GArray *q, vec_t *aux)
{
    char d6buf[MAXLINE];

    matrix_to_d6_canon(aux, dim, d6buf);
    key128_t k; d6_to_key128(d6buf, &k);

    if ( g_bucket_lookup(b, &k)!=NULL ) {
        return;
    }
    if (g_bucket_insert_copy128(b, &k)) {
        g_array_append_val(q, k);
    }
}

/**
 * populate_orbit - Populates the orbit of a given code by exploring its transformations.
 * @bucket: Pointer to a GHashBucket used for storing unique codes.
 * @code: The initial code to start orbit generation (as a string).
 * @mat: Pointer to a matrix vector used for transformations.
 * @aux: Auxiliary vector used during transformations.
 *
 * This function generates the orbit of the given code by applying column and row
 * transformations to its matrix representation. It uses a queue to explore all reachable
 * codes, adding new codes to the orbit if they result from valid transformations.
 * The function does not add the seed code to the bucket, only its transformations.
 *
 * Returns: The original code pointer if successful, or NULL if the input code is NULL.
 */
static char* populate_orbit(GHashBucket *bucket, char *code, vec_t *mat, vec_t *aux)
{
    if (code == NULL) {
        return NULL;
    }

    key128_t seedk; d6_to_key128(code, &seedk);

    GArray *q = g_array_new(FALSE, FALSE, sizeof(key128_t));
    guint   h = 0;

    g_array_append_val(q, seedk);
    // Do not add seedk to the bucket - no need

    char d6buf[MAXLINE];

    while (h < q->len) {
        key128_t cur = g_array_index(q, key128_t, h++);
        key128_to_d6(&cur, d6buf);

        matrix_from_d6(d6buf, mat, dim);

        for (ind_t i=0; i<dim; ++i) {
            conditional_add_col(mat, aux, dim, i);
            add_code(bucket, q, aux);
        }
        for (ind_t i=0; i<dim; ++i) {
            for (ind_t j=i+1; j<dim; ++j) {
                if (conditional_add_row(mat, aux, dim, i, j)) {
                    add_code(bucket, q, aux);
                }
                if (conditional_add_row(mat, aux, dim, j, i)) {
                    add_code(bucket, q, aux);
                }
            }
        }
    }
    printlog(3, "%s queue length: %u", code, q->len);
    g_array_free( q, TRUE );

    return code;
}

typedef struct {
    GHashBucket *bucket;
    gboolean    run;
    gulong      time;
    gulong      lines;
    gulong      reps;
} GBucketThreadData;

gpointer print_code_set_size(gpointer thread_data)
{
    GBucketThreadData *data = (GBucketThreadData*)thread_data;
    while (data->run) {
        printlog(2, "lines: %.2fM; reps: %.2fM; bucket: %.2fM", data->lines/1000000.0, data->reps/1000000.0, g_bucket_size(data->bucket)/1000000.0);
        g_usleep(data->time);
    }
    return NULL;
}

int main(int argc, char *argv[])
{
    tic();

    int opt = 1, v = 0;

    gulong t = 1;
    gsize  n = 1023;
    gsize  m = 0;

    FILE *in = stdin, *out = stdout;

    while ((opt = getopt(argc, argv, "vhn:i:o:t:m:")) != -1) {
        switch (opt) {
        case 'v':
            ++v;
            break;
        case 'n':
            n = (ind_t)atoi(optarg);
            break;
        case 't':
            t = (gulong)atoi(optarg);
            break;
        case 'm':
            m = (gsize)atoi(optarg);
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
        case 'h':
            help(argv[0]);
            exit(EXIT_SUCCESS);
        default: /* '?' */
            help(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    verbosity_level = v;

    char line[MAXLINE];
    if ( fgets(line, MAXLINE, in)==NULL ) {
        fprintf(stderr, "error reading input\n");
        // g_bucket_destroy(code_set);
        exit(1);
    }

    dim = graphsize( line );
    if (dim == 0) {
        fprintf(stderr, "read dimension zero, quitting...\n");
        // g_bucket_destroy(code_set);
        exit(1);
    }
    init_nauty_data( dim );

    vec_t *mat = init( dim );
    if (mat == NULL) {
        fprintf(stderr, "mat allocation error\n");
        // g_bucket_destroy( code_set );
        exit(1);
    }
    vec_t *aux = init( dim );
    if (aux == NULL) {
        fprintf(stderr, "mat allocation error\n");
        // g_bucket_destroy( code_set );
        exit(1);
    }

    char d6[MAXLINE];
    // digraph6_to_matrix(line, mat, dim);
    d6_to_d6_canon(line, d6);

    /* this auto shard calculation seems not to work as expected 
    if (m != 0) {
        n = shards_for_cap_limit( m*1000000ULL, 1.15, 0.80, 21 );
        printlog(1, "peak value (-m) set, calculated number of shards: %u", n);
    }
    */
    // create a hash table bucket
    GHashBucket* code_set = g_bucket_new_128(
        g_free,        // Function to free the key when the bucket is destroyed
        NULL,          // Function to free the value (not needed for simple booleans)
        n
    );
    g_bucket_reserve(code_set, m*1000000ULL);

    populate_orbit(code_set, d6, mat, aux);
    fprintf(out, "%s\n", d6);

    //gsize reps = 1;

    GBucketThreadData thread_data = { code_set, TRUE, t*1000000, 1, 1 };
    GThread *print_thread;
    if (v>=2) {
        print_thread = g_thread_new(NULL, print_code_set_size, &thread_data);
    }

    while (fgets(line, MAXLINE, in) != NULL) {
        ++thread_data.lines;
        // transform to canonical form
        d6_to_d6_canon(line, d6);
        key128_t repkey; d6_to_key128(d6, &repkey);
        // try to delete; if succesful, then continue
        if ( g_bucket_remove(code_set, &repkey) ) {
            continue;
        }
        // save the code
        // line[strlen(line)-1] = 0;
        fprintf(out, "%s\n", d6);
        ++thread_data.reps;
        // populate orbit
        populate_orbit(code_set, d6, mat, aux);
    }

    printlog(1, "%u representatives found", thread_data.reps);

    if (v>=2) {
        thread_data.run = FALSE;
        g_thread_join( print_thread );
    }

    free_nauty_data();
    free( mat );
    free( aux );

    // final cleaning up
    if (code_set) {
        g_bucket_destroy(code_set);
    }
    if (in != stdin) {
        fclose(in);
    }
    if (out != stdout) {
        fclose(out);
    }
    return 0;
}

