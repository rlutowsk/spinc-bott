#include "common.h"
#include "adjpack11.h"
#include "bucket.h"

static ADJPACK_INLINE void adjpack_from_matrix_legacy(const vec_t *mat, unsigned n, key128_t *out) {
    assert(n >= 1 && n <= 11);
#if D6PACK_HAVE_UINT128
	out->u = 0;
    uint64_t row;
    for (unsigned i = 0; i < n; ++i) {
        row = bitreverse64(mat[i]) >> (64 - n);
        out->u <<= n;
        out->u |= row;
    }
#else
    uint64_t lo = 0, hi = 0;
    unsigned pos = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (int j = n - 1; j >= 0; --j) {
            uint64_t bit = (mat[i] >> j) & 1;
            if (pos < 64)
                lo |= bit << pos;
            else
                hi |= bit << (pos - 64);
            ++pos;
        }
    }
    memcpy(out->b + 0, &lo, 8);
    memcpy(out->b + 8, &hi, 8);
#endif
    adjpack_zero_above_nsquare(out, n);
    adjpack_set_n(out, n);
}

static ADJPACK_INLINE void adjpack_from_matrix_native(const vec_t *mat, unsigned n, key128_t *out) {
    assert(n >= 1 && n <= 11);
#if D6PACK_HAVE_UINT128
    // __uint128_t row = 0;
    out->u = 0;
    unsigned j = 0;
    for (unsigned i = 0; i < n; ++i, j += n) {
        out->u |= (__uint128_t)mat[i] << (j);
    }
    out->u = bitreverse128(out->u) >> (128 - j);
#else
    uint64_t lo = 0, hi = 0;
    unsigned pos = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (int j = n - 1; j >= 0; --j) {
            uint64_t bit = (mat[i] >> j) & 1;
            if (pos < 64)
                lo |= bit << pos;
            else
                hi |= bit << (pos - 64);
            ++pos;
        }
    }
    memcpy(out->b + 0, &lo, 8);
    memcpy(out->b + 8, &hi, 8);
#endif
    adjpack_zero_above_nsquare(out, n);
    adjpack_set_n(out, n);
}

bool test_correctness() {
    // printf("Testing correctness of adjpack_from_matrix vs adjpack_from_matrix_legacy...\n");

    for (unsigned dim = 1; dim <= 11; ++dim) {
        vec_t *mat = init(dim);
        bool all_passed = true;

        // Test dla kilku przykładowych macierzy
        for (unsigned i = 0; i < 100; ++i) {
            for (unsigned r = 0; r < dim; ++r) {
                mat[r] = i ^ (r * 0xA5A5A5A5A5A5A5A5ULL); // przykładowe dane
                mat[r] &= ((1ULL << dim) - 1); // ograniczenie do dim bitów
            }

            key128_t a, b;
            adjpack_from_matrix_native(mat, dim, &a);
            adjpack_from_matrix_legacy(mat, dim, &b);

            if (!key128_equal(&a, &b)) {
                printf("    Mismatch at dim=%u, test=%u\n", dim, i);
                all_passed = false;
                break;
            }
        }

        if (all_passed)
            printf("    dim=%u: OK\n", dim);
        else
            return false;
    }
    return true;
}


double get_time_ms() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}


void benchmark(bool verbose) {
    const unsigned dim = 10;
    const unsigned N = 0xFFFFFFF0; //1<<31;
    vec_t *mat = init(dim);
    volatile key128_t out;

    // Wypełnij przykładową macierz
    for (unsigned r = 0; r < dim; ++r)
        mat[r] = (r * 0x123456789ABCDEFULL) & ((1ULL << dim) - 1);

    double start, end;

    start = get_time_ms();
    for (unsigned i = 0; i < N; ++i)
        adjpack_from_matrix_native(mat, dim, (key128_t*)&out);
    end = get_time_ms();
    double t1 = end - start;

    start = get_time_ms();
    for (unsigned i = 0; i < N; ++i)
        adjpack_from_matrix_legacy(mat, dim, (key128_t*)&out);
    end = get_time_ms();
    double t2 = end - start;

    printf("    Benchmark for dim=%u, %u iterations:\n", dim, N);
    printf("      adjpack_from_matrix:         %.6f ms\n", t1);
    printf("      adjpack_from_matrix_legacy:  %.6f ms\n", t2);

    if (verbose) {
        printf("      Speed ratio (legacy / 128): %.2f%%\n", (t2 / t1) * 100.0);
    }
}

int main() {
    printf("=== [adjpack_from_matrix] testing implementations ===\n");
    fflush(stdout);
    if (!test_correctness()) {
        printf("=== [adjpack_from_matrix] some tests failed ===\n");
        return 1;
    }
    printf("=== [adjpack_from_matrix] all tests passed ===\n"); fflush(stdout);
    benchmark(true);
    printf("=== [adjpack_from_matrix] benchmark completed ===\n");
    return 0;
}
