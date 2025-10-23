#include "common.h"
#include "adjpack11.h"
#include <inttypes.h>

#if defined(__clang__) && __clang_major__ >= 5
static ADJPACK_INLINE uint64_t bitreverse64_legacy(uint64_t x) {
    // Check if the compiler supports a built-in function for bit reversal.
    // Clang supports this from version 5.
    return __builtin_bitreverse64(x);
}
#else
static ADJPACK_INLINE uint64_t bitreverse64_legacy(uint64_t x) {
    // Portable, and very efficient implementation as a fallback.
    x = ((x & 0x5555555555555555ULL) <<  1) | ((x >>  1) & 0x5555555555555555ULL);
    x = ((x & 0x3333333333333333ULL) <<  2) | ((x >>  2) & 0x3333333333333333ULL);
    x = ((x & 0x0F0F0F0F0F0F0F0FULL) <<  4) | ((x >>  4) & 0x0F0F0F0F0F0F0F0FULL);
    x = ((x & 0x00FF00FF00FF00FFULL) <<  8) | ((x >>  8) & 0x00FF00FF00FF00FFULL);
    x = ((x & 0x0000FFFF0000FFFFULL) << 16) | ((x >> 16) & 0x0000FFFF0000FFFFULL);
    x = (x << 32) | (x >> 32);
    return x;
}
#endif

static ADJPACK_INLINE __uint128_t bitreverse128_native(__uint128_t x) {
#if defined(__clang__) && __clang_major__ >= 18
    return __builtin_bitreverse128(x);
#else
    // Maski 128-bitowe zbudowane z dwóch 64-bitowych części
    const __uint128_t mask_1 = ((__uint128_t)0x5555555555555555ULL << 64) | 0x5555555555555555ULL;
    const __uint128_t mask_2 = ((__uint128_t)0x3333333333333333ULL << 64) | 0x3333333333333333ULL;
    const __uint128_t mask_4 = ((__uint128_t)0x0F0F0F0F0F0F0F0FULL << 64) | 0x0F0F0F0F0F0F0F0FULL;

    // Trzy kroki odwracania bitów
    x = ((x & mask_1) << 1) | ((x >> 1) & mask_1);
    x = ((x & mask_2) << 2) | ((x >> 2) & mask_2);
    x = ((x & mask_4) << 4) | ((x >> 4) & mask_4);

    // Odwrócenie bajtów
    #if defined(__GNUC__) && __GNUC__ >= 12
        return __builtin_bswap128(x);
    #else
        uint64_t hi = (uint64_t)(x >> 64);
        uint64_t lo = (uint64_t)x;
        hi = __builtin_bswap64(hi);
        lo = __builtin_bswap64(lo);
        return ((__uint128_t)lo << 64) | hi;
    #endif
#endif
}

static ADJPACK_INLINE __uint128_t bitreverse128_legacy(__uint128_t x) {
    // Portable, and very efficient implementation as a fallback.
    uint64_t hi = (uint64_t)(x >> 64);
    uint64_t lo = (uint64_t)x;
    
    return bitreverse64(hi) | ((__uint128_t)bitreverse64(lo) << 64);
}

bool test_correctness() {
    uint64_t test_values[] = {
        0x0000000000000001ULL,
        0x8000000000000000ULL,
        0x0123456789ABCDEFULL,
        0xFFFFFFFFFFFFFFFFULL,
        0x00000000FFFFFFFFULL,
        0xF0F0F0F0F0F0F0F0ULL
    };

    for (int i = 0; i < sizeof(test_values)/sizeof(test_values[0]); ++i) {
        uint64_t x = test_values[i];
        uint64_t r1 = bitreverse64(x);
        uint64_t r2 = bitreverse64_legacy(x);
        if (r1 != r2) {
            printf("    Test %d FAILED: input = 0x%016" PRIx64 "\n", i, x);
            printf("      bitreverse64        = 0x%016" PRIx64 "\n", r1);
            printf("      bitreverse64_legacy = 0x%016" PRIx64 "\n", r2);
            return false;
        }
    }
    printf("    bitreverse64:  OK\n");
    return true;
}

void benchmark(bool show) {
    const uint64_t N = 1UL << 24;
    volatile uint64_t sum = 0;
    clock_t start, end;

    start = clock();
    for (uint64_t i = 0; i < N; ++i)
        sum += bitreverse64(i);
    end = clock();
    if( show )
        printf("    bitreverse64:         %6.2f ms\n", 1000.0 * (end - start) / CLOCKS_PER_SEC);

    sum = 0;
    start = clock();
    for (uint64_t i = 0; i < N; ++i)
        sum += bitreverse64_legacy(i);
    end = clock();
    if( show )
        printf("    bitreverse64_legacy:  %6.2f ms\n", 1000.0 * (end - start) / CLOCKS_PER_SEC);
}

bool test_bitreverse128_correctness() {
    __uint128_t tests[] = {
        (__uint128_t)1,
        (__uint128_t)1 << 127,
        (__uint128_t)0xFFFFFFFFFFFFFFFFULL,
        (__uint128_t)0xFFFFFFFFFFFFFFFFULL << 64,
        ((__uint128_t)0x0123456789ABCDEFULL << 64) | 0xFEDCBA9876543210ULL,
        ((__uint128_t)0x195D3B7F084C2A6EULL << 64) | 0xF7B3D591E6A2C480ULL,
        (__uint128_t)0,
        (__uint128_t)-1,
    };

    for (int i = 0; i < sizeof(tests)/sizeof(tests[0]); ++i) {
        __uint128_t result = bitreverse128_native(tests[i]);
        __uint128_t expected = bitreverse128_legacy(tests[i]);
        if (result != expected) {
            printf("    bitreverse128 test %d FAILED\n", i);
            printf("      input    = 0x%016lx%016lx\n", (uint64_t)(tests[i] >> 64), (uint64_t)tests[i]);
            printf("      expected = 0x%016lx%016lx\n", (uint64_t)(expected >> 64), (uint64_t)expected);
            printf("      result   = 0x%016lx%016lx\n", (uint64_t)(result >> 64), (uint64_t)result);
            return false;
        }
    }
    printf("    bitreverse128: OK\n");
    return true;
}

void benchmark128() {
    const int N = 1<<26;
    volatile __uint128_t sum = 0;
    clock_t start, end;

    // Rozgrzewka
    for (int i = 0; i < N; ++i) {
        __uint128_t x = ((__uint128_t)i << 64) | (N - i);
        sum += bitreverse128_native(x);
    }

    // Benchmark bitreverse128
    sum = 0;
    start = clock();
    for (int i = 0; i < N; ++i) {
        __uint128_t x = ((__uint128_t)i << 64) | (N - i);
        sum += bitreverse128_native(x);
    }
    end = clock();
    printf("    bitreverse128:        %6.2f ms\n", 1000.0 * (end - start) / CLOCKS_PER_SEC);

    // Benchmark bitreverse128_legacy
    sum = 0;
    start = clock();
    for (int i = 0; i < N; ++i) {
        __uint128_t x = ((__uint128_t)i << 64) | (N - i);
        sum += bitreverse128_legacy(x);
    }
    end = clock();
    printf("    bitreverse128_legacy: %6.2f ms\n", 1000.0 * (end - start) / CLOCKS_PER_SEC);
}

int main() {
    printf("=== [bitreverse] testing 64/128 implementations ===\n");
    fflush(stdout);
    if (!test_correctness() || !test_bitreverse128_correctness()) {
        printf("=== [bitreverse] some tests failed ===\n");
        return 1;
    }
    printf("=== [bitreverse] all tests passed ===\n");
    fflush(stdout);
    benchmark(false);
    benchmark(true);

    benchmark128();
    printf("=== [bitreverse] benchmark completed ===\n");
    return 0;
}
