#include "bott.h"
#include "common.h"
#include "dag.h"    /* matrix_to_d6 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <time.h>
#include <errno.h>

/* --- Pomoc dla CLI --- */
static void help(const char *name)
{
    fprintf(stderr,
        "Usage: %s [-j njobs] [-s start_dim] [-d dimension] [-a] [-p] [-n] [-v] [-h] [-o <path>|stdout|-]\n"
        "-j: number of threads\n"
        "-d: target dimension to calculate\n"
        "-s: starting dimension, between 3 and 11, less than target dimension\n"
        "-a: calculate spin structures also\n"
        "-p: show progress during computation\n"
        "-n: suppress final numeric output\n"
        "-o: write DAG codes (d6) to file; use '-' or 'stdout' to write to STDOUT\n"
        "-v: be verbose\n",
        name);
}

/* --- Z oryginału: get_max_state() --- */
static inline state_t get_max_state(ind_t dim)
{
    ind_t c;
    state_t max_state;
    for (c = 0, max_state = 1; c < dim - 1; c++) {
        max_state <<= (dim - c - 2);
    }
    return max_state - 1;
}

/* global */
int calculate_spin = 0;

/* --- I/O bufor per‑task dla kodów d6 --- */
struct out_ctx {
    int enabled;        /* czy włączono -o */
    FILE *fp;           /* dokąd pisać (stdout/plik) */
    char *buf;          /* bufor dla batchowego I/O */
    size_t used;        /* ile bajtów zajęte */
    size_t cap;         /* pojemność bufora */
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

static inline void increase_dimension(vec_t *mat, ind_t dim)
{
    for (ind_t i = 1; i < dim; i++) { mat[i] <<= 1; }
}
static inline void decrease_dimension(vec_t *mat, ind_t dim)
{
    for (ind_t i = 1; i < dim; i++) { mat[i] >>= 1; }
}

/* --- Deklaracje --- */
size_t backtrack(vec_t *mat, vec_t **cache, ind_t cdim, ind_t ddim,
                 size_t *spinc, size_t *spin, struct out_ctx *out);

/* --- Monitor postępu: drukuje na wskazany strumień (stderr przy -o) --- */
static void monitor_progress(unsigned long *progress_ptr, unsigned long total, FILE *stream)
{
    const double period = 0.25;
    double next = omp_get_wtime();
    struct timespec ts = { .tv_sec = 0, .tv_nsec = 50 * 1000 * 1000 }; /* 50 ms */

    for (;;) {
        unsigned long cur;
        #pragma omp atomic read
        cur = *progress_ptr;
        if (cur >= total) break;

        double now = omp_get_wtime();
        if (now >= next) {
            double pct = (total == 0) ? 100.0 : (100.0 * (double)cur / (double)total);
            fprintf(stream, "\33[2K\r%8.4f%%", pct);
            fflush(stream);
            next = now + period;
        }
        nanosleep(&ts, NULL);
    }
}

/* --- Właściwe backtrack: identyczna logika jak w Twojej wersji, + opcjonalny zapis kodów d6 --- */
size_t backtrack(vec_t *mat, vec_t **cache, ind_t cdim, ind_t ddim,
                 size_t *spinc, size_t *spin, struct out_ctx *out)
{
    vec_t max = 1 << (cdim - 2);
    vec_t r;
    vec_t row = ddim - cdim;

    increase_dimension(mat, ddim);

    if (row == 0) {
        /*
         * final step:
         * nie trzeba liczyć macierzy [0 000; 0 mat]; spinc info z mat identyczne
         */
        char code_buf[256];

        mat[0] = 0;
        *spinc += 1;
        if (out && out->enabled) {
            matrix_to_d6(mat, ddim, code_buf);
            out_append_line(out, code_buf);
        }

        if (calculate_spin) {
            *spin += is_spin(mat, ddim);
            for (r = 1; r < max; r++) {
                mat[0] = cache[ddim][r];
                if (is_spinc(mat, ddim)) {
                    *spinc += 1;
                    *spin  += is_spin(mat, ddim);
                    if (out && out->enabled) {
                        matrix_to_d6(mat, ddim, code_buf);
                        out_append_line(out, code_buf);
                    }
                }
            }
        } else {
            for (r = 1; r < max; r++) {
                mat[0] = cache[ddim][r];
                if (is_spinc(mat, ddim)) {
                    *spinc += 1;
                    if (out && out->enabled) {
                        matrix_to_d6(mat, ddim, code_buf);
                        out_append_line(out, code_buf);
                    }
                }
            }
        }
    } else {
        /* rekurencja */
        for (r = 0; r < max; r++) {
            mat[row] = cache[cdim][r];
            if (is_spinc(&mat[row], cdim)) {
                backtrack(mat, cache, cdim + 1, ddim, spinc, spin, out);
                /* po każdym zejściu zmniejszamy wymiar, bo następne wołanie znów go zwiększy */
                decrease_dimension(mat, ddim);
            }
        }
    }
    return *spinc;
}

/* -------------------- main -------------------- */
int main(int argc, char *argv[])
{
    int opt = 1, bs = -1, no_output = 0;
    state_t test = 0;
    int progress_enabled = 0;

    /* domyślne wartości */
    ind_t sdim = 0, dim = 6;
    vec_t *cache[12] = { 0,0,0,0,0,0,0,0,0,0,0,0 };

    state_t loop_start = 1, loop_stop = 1;
    vec_t row;
    size_t spinc = 0, spin = 0;
    size_t cache_size;

    /* -o */
    int output_enabled = 0;
    const char *out_path = NULL;
    FILE *out_fp = NULL;

    tic();

    while ((opt = getopt(argc, argv, "vhj:d:s:ac:npt:o:")) != -1) {
        switch (opt) {
        case 't': test = atol(optarg); break;
        case 'p': progress_enabled = 1; break;
        case 'n': no_output = 1; break;
        case 'c': bs = atoi(optarg); break;
        case 'a': calculate_spin = 1; break;
        case 'v': verbosity_level++; break;
        case 'j': omp_set_num_threads(atoi(optarg)); break;
        case 'd': dim  = (ind_t)atoi(optarg); break;
        case 's': sdim = (ind_t)atoi(optarg); break;
        case 'o': output_enabled = 1; out_path = optarg; break;
        case 'h': help(argv[0]); exit(EXIT_SUCCESS);
        default: help(argv[0]); exit(EXIT_FAILURE);
        }
    }

    if (sdim == 0) {
        sdim = (dim > 11) ? 11 : dim - 1;
    }
    if (sdim < 3 || sdim > 11 || sdim >= dim) {
        fprintf(stderr, "Starting dimension (%i) out of range, dim is %i.\n", sdim, dim);
        exit(EXIT_FAILURE);
    }

    /* cache */
    for (int i = sdim; i < 12; i++) {
        vec_t *tmp;
        populate_cache(&tmp, &cache_size, i);
        cache[i] = tmp;
    }

    state_t max_state = get_max_state(sdim);
    row = dim - sdim;

    /* heurystyka chunk size (historyczna – nie używana bezpośrednio w tasks, ale zostawiona) */
    if (bs == -1) {
        bs = 7 * sdim - 49;
        bs = (bs < 0) ? 0 : bs;
    }

    if (test > 0) {
        loop_start = (max_state + 1) / 2;
        loop_stop  = loop_start + ((state_t)1 << test) - 1;
    } else {
        loop_start = 0;
        loop_stop  = max_state;
    }

    /* -o: otwarcie wyjścia */
    if (output_enabled) {
        if (strcmp(out_path, "-") == 0 || strcasecmp(out_path, "stdout") == 0) {
            out_fp = stdout;
        } else {
            out_fp = fopen(out_path, "w");
            if (!out_fp) {
                fprintf(stderr, "Cannot open output file '%s': %s\n", out_path, strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
    }

    const unsigned long total_iters = (unsigned long)(loop_stop - loop_start + 1);

    /* dobór ziarnistości tasków: ~8 zadań na wątek, z bezpiecznymi limitami */
    int nthreads = omp_get_max_threads();
    unsigned long target_tasks = (unsigned long)(nthreads > 0 ? nthreads : 1) * 8ul;
    unsigned long grain = (total_iters + target_tasks - 1) / (target_tasks ? target_tasks : 1);
    if (grain < 256ul)   grain = 256ul;
    if (grain > 65536ul) grain = 65536ul;

    printlog(2, "iteration range: %lu..%lu (total=%lu)\n",
             (unsigned long)loop_start, (unsigned long)loop_stop, total_iters);
    printlog(2, "openmp task grain set to %lu iterations\n", grain);

    /* postęp i strumienie */
    static unsigned long progress = 0;
    #pragma omp atomic write
    progress = 0;

    FILE *progress_stream = output_enabled ? stderr : stdout;
    FILE *summary_stream  = (output_enabled ? stderr : stdout);

    /* ------------------ Równoległość (zadania) ------------------ */
    #pragma omp parallel
    #pragma omp single nowait
    {
        /* monitor postępu */
        if (progress_enabled) {
            #pragma omp task untied firstprivate(total_iters, progress_stream)
            {
                monitor_progress(&progress, total_iters, progress_stream);
            }
        }

        /* zadania robocze */
        #pragma omp taskgroup
        {
            for (unsigned long base = (unsigned long)loop_start; base <= (unsigned long)loop_stop; ) {
                unsigned long end = base + grain - 1;
                if (end > (unsigned long)loop_stop) end = (unsigned long)loop_stop;

                #pragma omp task firstprivate(base, end, row, dim, sdim) shared(cache, calculate_spin, progress, out_fp, output_enabled) \
                                 untied
                {
                    vec_t *tmat = init(dim);
                    size_t local_spinc = 0, local_spin = 0;

                    struct out_ctx out = {0};
                    out.enabled = output_enabled;
                    out.fp = out_fp ? out_fp : stdout;
                    if (out.enabled) {
                        out.cap = OUTBUF_CAP;
                        out.buf = (char*)malloc(out.cap);
                        out.used = 0;
                        if (!out.buf) {
                            fprintf(stderr, "malloc failed for output buffer\n");
                            exit(EXIT_FAILURE);
                        }
                    }

                    for (state_t s = (state_t)base; s <= (state_t)end; ++s) {
                        matrix_by_state(&tmat[row], cache[sdim], s, sdim);
                        if (is_spinc(&tmat[row], sdim)) {
                            backtrack(tmat, cache, sdim + 1, dim, &local_spinc, &local_spin, out.enabled ? &out : NULL);
                            decrease_dimension(tmat, dim);
                        }
                    }

                    /* spłucz bufor kodów i sprzątnij */
                    if (out.enabled) {
                        out_flush(&out);
                        free(out.buf);
                    }
                    free(tmat);

                    /* aktualizacja liczników globalnych */
                    #pragma omp atomic update
                    spinc += local_spinc;
                    if (calculate_spin) {
                        #pragma omp atomic update
                        spin += local_spin;
                    }

                    /* postęp: cały chunk jednym atomikiem */
                    unsigned long done = (unsigned long)(end - base + 1);
                    #pragma omp atomic update
                    progress += done;
                }

                if (end == (unsigned long)loop_stop) break;
                base = end + 1;
            }
        }

        /* domknij postęp */
        #pragma omp atomic write
        progress = total_iters;
    }

    if (progress_enabled) {
        fprintf(progress_stream, "\33[2K\r%8.4f%%\n", 100.0);
        fflush(progress_stream);
    }

    /* Log końcowy */
    if (calculate_spin) {
        printlog(1, "calculations finished: %lu/%lu spin/spinc manifolds in total\n",
                 (unsigned long)spin, (unsigned long)spinc);
    } else {
        printlog(1, "calculations finished: %lu spinc manifolds in total\n",
                 (unsigned long)spinc);
    }

    /* Finalne podsumowanie liczbowe (stderr przy -o, aby nie mieszać się z kodami na stdout/plik) */
    if (!no_output) {
        if (calculate_spin) {
            fprintf(summary_stream, "%lu/%lu\n", (unsigned long)spin, (unsigned long)spinc);
        } else {
            fprintf(summary_stream, "%lu\n", (unsigned long)spinc);
        }
        fflush(summary_stream);
    }

    /* zamknij plik, jeśli to nie stdout */
    if (output_enabled && out_fp && out_fp != stdout) {
        fclose(out_fp);
    }

    /* sprzątanie cache */
    for (int i = sdim; i < 12; i++) { free(cache[i]); }

    return 0;
}
