#include "common.h"

#include "bott.h"
#include "dag.h"    /* matrix_to_d6 */
#include "adjpack11.h"
#include "bucket.h"
#include "tlsbuf.h"

#include <stdlib.h>
#include <string.h>
// #include <unistd.h>
#include <omp.h>
// #include <errno.h>

/* --- Pomoc dla CLI --- */
static void help(const char *name)
{
    fprintf(stderr,
        "Usage: %s [-j njobs] [-s start_dim] [-d dimension] [-a] [-p] [-n] [-v] [-h] [-o <path>|stdout|-]\n"
        "-j: number of threads\n"
        "-d: dimension to calculate\n"
        "-p: show progress during computation\n"
        "-o: write DAG codes (d6) to file; stdout is default\n"
        "-v: be verbose\n",
        name);
}

/* --- I/O buffer per‑task for d6 codes --- */
struct out_ctx {
    int enabled;        /* whether -o is enabled */
    FILE *fp;           /* where to write (stdout/file) */
    char *buf;          /* buffer for batched I/O */
    size_t used;        /* how many bytes are used */
    size_t cap;         /* buffer capacity */
};
#define OUTBUF_CAP (1u<<20) /* ~1 MiB per task */

static inline void out_flush(struct out_ctx *out)
{
    if (!out || !out->enabled || out->used == 0) return;
    #pragma omp critical(output)
    {
        (void)fwrite(out->buf, 1, out->used, out->fp);
        fflush(out->fp);
    }
    out->used = 0;
}

static inline void out_append_line(struct out_ctx *out, const char *line)
{
    if (!out || !out->enabled) return;
    size_t len = strlen(line);
    if (out->used + len + 1 >= out->cap) {
        out_flush(out);
    }
    memcpy(out->buf + out->used, line, len);
    out->used += len;
    out->buf[out->used++] = '\n';
}

/* -------------------- main -------------------- */
int main(int argc, char *argv[])
{
    int opt = 1;
    int progress_enabled = 0;
    /* defaults */
    ind_t dim = 6;
    vec_t *cache;
    size_t cache_size;

    /* -o */
    const char *out_path = NULL;
    FILE *out_fp = stdout;

    tic();

    while ((opt = getopt(argc, argv, "vhj:d:po:")) != -1) {
        switch (opt) {
        case 'p': progress_enabled = 1; break;
        case 'v': verbosity_level++; break;
        case 'j': omp_set_num_threads(atoi(optarg)); break;
        case 'd': dim  = (ind_t)atoi(optarg); break;
        case 'o': out_path = optarg; break;
        case 'h': help(argv[0]); exit(EXIT_SUCCESS);
        default: help(argv[0]); exit(EXIT_FAILURE);
        }
    }

    /* cache */
    populate_cache(&cache, &cache_size, dim);

    state_t max_state = get_max_state(dim);

    /* -o: open output */
    if (out_path!=NULL) {
        out_fp = fopen(out_path, "w");
        if (!out_fp) {
            fprintf(stderr, "Cannot open output file '%s': %s\n", out_path, strerror(errno));
            exit(EXIT_FAILURE);
        }
    }

    const unsigned long total_iters = (unsigned long)(max_state + 1);

    /* task granularity selection: ~8 tasks per thread, with safe limits */
    int nthreads = omp_get_num_threads();
    unsigned long target_tasks = (unsigned long)(nthreads > 0 ? nthreads : 1) * 8ul;
    unsigned long grain = (total_iters + target_tasks - 1) / (target_tasks ? target_tasks : 1);
    if (grain < 256ul)   grain = 256ul;
    if (grain > 131072ul) grain = 131072ul;

    printlog(2, "iteration range: %lu..%lu (total=%lu)\n",
             0lu, (unsigned long)max_state, total_iters);
    printlog(2, "openmp task grain set to %lu iterations\n", grain);

    /* progress and streams */
    static unsigned long progress = 0;
    #pragma omp atomic write
    progress = 0;

    FILE *progress_stream = stderr;

    GHashBucket *g_canonical_set = g_bucket_new_128(NULL, NULL, 1023);

    printlog(1, "starting calculations");
    /* ------------------ Parallelism (tasks) ------------------ */
    #pragma omp parallel
    #pragma omp single nowait
    {
        /* progress monitor */
        if (progress_enabled) {
            #pragma omp task untied firstprivate(total_iters, progress_stream)
            {
                monitor_progress(&progress, total_iters, progress_stream);
            }
        }

        /* worker tasks */
        #pragma omp taskgroup
        {
            for (unsigned long base = 0ul; base <= (unsigned long)max_state; ) {
                unsigned long end = base + grain - 1;
                if (end > (unsigned long)max_state) end = (unsigned long)max_state;

                #pragma omp task firstprivate(base, end, dim) shared(cache, progress, out_fp) \
                                 untied
                {
                    vec_t *mat = init(dim),
                          *can = init(dim);
                    key128_t key;
                    char buf[MAXLINE];

                    struct out_ctx out = {0};
                    out.enabled = true;
                    out.fp = out_fp;
                    out.cap = OUTBUF_CAP;
                    out.buf = (char*)malloc(out.cap);
                    out.used = 0;
                    if (!out.buf) {
                        fprintf(stderr, "malloc failed for output buffer\n");
                        exit(EXIT_FAILURE);
                    }
                    
                    init_nauty_data(dim);
                    for (state_t s = (state_t)base; s <= (state_t)end; ++s) {
                        matrix_by_state(mat, cache, s, dim);
                        assert(is_orientable(mat, dim));
                        // matrix_to_matrix_canon(mat, dim, can);
                        // adjpack_from_matrix(can, dim, &key);
                        matrix_to_key128_canon(mat, dim, &key);
                        if (g_bucket_insert_copy128(g_canonical_set, &key)) {
                            d6pack_encode(&key, dim, buf);
                            out_append_line(&out, buf);
                        }
                    }
                    free_nauty_data();

                    /* flush the code buffer and clean up */
                    out_flush(&out);
                    free(out.buf);

                    /* progress: entire chunk with a single atomic update */
                    unsigned long done = (unsigned long)(end - base + 1);
                    #pragma omp atomic update
                        progress += done;
                    
                    free(mat);
                    free(can);
                }

                if (end == (unsigned long)max_state) break;
                base = end + 1;
            }
        }

        /* finalize progress */
        #pragma omp atomic write
            progress = total_iters;
    }

    if (progress_enabled) {
        fprintf(progress_stream, "\33[2K\r%8.4f%%\n", 100.0);
        fflush(progress_stream);
    }
    printlog(1, "all calculations done");

    g_bucket_destroy(g_canonical_set);

    /* close the file if it's not stdout */
    if (out_fp != stdout) {
        fclose(out_fp);
    }

    /* clean up cache */
    free(cache);

    return 0;
}
