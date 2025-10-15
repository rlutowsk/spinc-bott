#include "bott.h"
#include "common.h"
#include "dag.h"
#include "bucket.h"
#include "adjpack11.h"

#define MAXLINE 128

#include <assert.h>

#define INITIAL_CAPACITY 1024

typedef struct {
    vec_t *rows;     // flat array of rows
    size_t len;      // number of matrices stored
    size_t capacity; // total number of matrices that can be stored
    size_t dim;      // number of rows per matrix
} MatArray;

// Create a new MatArray
MatArray *matarray_create(size_t dim) {
    MatArray *arr = malloc(sizeof(MatArray));
    arr->dim = dim;
    arr->len = 0;
    arr->capacity = INITIAL_CAPACITY;
    arr->rows = malloc(sizeof(vec_t) * dim * arr->capacity);
    return arr;
}

// Free the MatArray
void matarray_free(MatArray *arr) {
    if (arr) {
        free(arr->rows);
        free(arr);
    }
}

// Append a matrix (vec_t[dim]) to the array
void matarray_append(MatArray *arr, const vec_t *mat) {
    if (arr->len >= arr->capacity) {
        arr->capacity *= 2;
        arr->rows = realloc(arr->rows, sizeof(vec_t) * arr->dim * arr->capacity);
        assert(arr->rows != NULL);
    }
    memcpy(&arr->rows[arr->len * arr->dim], mat, sizeof(vec_t) * arr->dim);
    arr->len++;
}

// Retrieve a pointer to the i-th matrix
vec_t *matarray_get(MatArray *arr, size_t index) {
    assert(index < arr->len);
    return &arr->rows[index * arr->dim];
}

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
static inline void add_code(GHashBucket *b, MatArray *q, vec_t *aux)
{
    /* 
    char d6buf[MAXLINE];

    matrix_to_d6_canon(aux, dim, d6buf);
    key128_t k; d6_to_key128(d6buf, &k);
    */
    // /*
    vec_t out[dim];
    matrix_to_matrix_canon(aux, dim, out);
    key128_t k; adjpack_from_matrix(out, dim, &k);
    // */
    if ( g_bucket_lookup(b, &k)!=NULL ) {
        return;
    }
    if (g_bucket_insert_copy128(b, &k)) {
        matarray_append(q, aux);
    }
}

/*
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
static char* populate_orbit(GHashBucket *bucket, char *code)
{
    if (code == NULL) {
        return NULL;
    }

    vec_t aux[dim];

    matrix_from_d6(code, aux, dim);

    MatArray *q = matarray_create(dim);

    matarray_append(q, aux);

    for (guint h=0; h < q->len; ++h) {
        vec_t *cur;

        for (ind_t i=0; i<dim; ++i) {
            cur = matarray_get(q, h);
            conditional_add_col(cur, aux, dim, i);
            add_code(bucket, q, aux);
        }
        for (ind_t i=0; i<dim; ++i) {
            for (ind_t j=i+1; j<dim; ++j) {
                cur = matarray_get(q, h);
                if (conditional_add_row(cur, aux, dim, i, j)) {
                    add_code(bucket, q, aux);
                }
                cur = matarray_get(q, h);
                if (conditional_add_row(cur, aux, dim, j, i)) {
                    add_code(bucket, q, aux);
                }
            }
        }
    }
    matarray_free( q );

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

    char d6[MAXLINE];
    // digraph6_to_matrix(line, mat, dim);
    d6_to_d6_canon(line, d6);

    // create a hash table bucket
    GHashBucket* code_set = g_bucket_new_128(
        g_free,        // Function to free the key when the bucket is destroyed
        NULL,          // Function to free the value (not needed for simple booleans)
        n
    );
    g_bucket_reserve(code_set, m*1000000ULL);

    populate_orbit(code_set, d6);
    fprintf(out, "%s\n", d6);

    GBucketThreadData thread_data = { code_set, TRUE, t*1000000, 1, 1 };
    GThread *print_thread;
    if (v>=2) {
        print_thread = g_thread_new(NULL, print_code_set_size, &thread_data);
    }

    while (fgets(line, MAXLINE, in) != NULL) {
        ++thread_data.lines;
        /* transform to canonical form */
        d6_to_d6_canon(line, d6);
        key128_t repkey;
        d6_to_key128(d6, &repkey);
        // try to delete; if succesful, then continue
        if ( g_bucket_remove(code_set, &repkey) ) {
            continue;
        }
        // save the code
        // line[strlen(line)-1] = 0;
        fprintf(out, "%s\n", d6);
        ++thread_data.reps;
        // populate orbit
        populate_orbit(code_set, d6);
    }

    printlog(1, "%u representatives found", thread_data.reps);

    if (v>=2) {
        thread_data.run = FALSE;
        g_thread_join( print_thread );
    }

    free_nauty_data();

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

