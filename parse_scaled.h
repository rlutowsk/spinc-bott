#pragma once

#include <ctype.h>

/* Parses decimal with optional suffix:
   - K/M/G/T/P/E   => *1000, *1000^2, ...
   - Ki/Mi/Gi/...  => *1024, *1024^2, ...
   Optional trailing 'B'/'b' is accepted (e.g., "2MiB").
   Returns true on success; false on invalid input or overflow. */
static INLINE bool parse_scaled_size(const char *s, size_t *out)
{
    if (!s || !*s || !out) return false;

    errno = 0;
    char *end = NULL;
    unsigned long long base = strtoull(s, &end, 10);
    if (errno != 0 || end == s) return false;   // no number parsed

    // Skip whitespace between number and suffix
    while (*end && isspace((unsigned char)*end)) ++end;

    unsigned long long mult = 1ULL;

    if (*end) {
        char c1 = (char)toupper((unsigned char)*end); 
        char c2 = (char)toupper((unsigned char)*(end + 1));
        const char *p = end;

        // Binary suffix: Ki, Mi, Gi, Ti, Pi, Ei (optional trailing 'B')
        if ((c1=='K'||c1=='M'||c1=='G'||c1=='T'||c1=='P'||c1=='E') && c2=='I') {
            unsigned power = 0;
            switch (c1) {
                case 'K': power = 1; break;
                case 'M': power = 2; break;
                case 'G': power = 3; break;
                case 'T': power = 4; break;
                case 'P': power = 5; break;
                case 'E': power = 6; break;
            }
            p += 2;
            if (*p=='B' || *p=='b') ++p;
            if (*p) return false;  // extra junk

            // mult = 1024^power, overflow‑safe
            mult = 1ULL;
            for (unsigned i = 0; i < power; ++i) {
                if (mult > ULLONG_MAX / 1024ULL) return false;
                mult *= 1024ULL;
            }
        }
        // Decimal suffix: K/M/G/T/P/E (optional trailing 'B')
        else if (c1=='K'||c1=='M'||c1=='G'||c1=='T'||c1=='P'||c1=='E') {
            switch (c1) {
                case 'K': mult = 1000ULL; break;
                case 'M': mult = 1000ULL * 1000ULL; break;
                case 'G': mult = 1000ULL * 1000ULL * 1000ULL; break;
                case 'T': mult = 1000ULL * 1000ULL * 1000ULL * 1000ULL; break;
                case 'P': mult = 1000ULL * 1000ULL * 1000ULL * 1000ULL * 1000ULL; break;
                case 'E': mult = 1000ULL * 1000ULL * 1000ULL * 1000ULL * 1000ULL * 1000ULL; break;
            }
            p += 1;
            if (*p=='B' || *p=='b') ++p;
            if (*p) return false;  // extra junk
        }
        else {
            // No recognized suffix: reject trailing junk
            if (*end) return false;
        }
    }

    // Overflow‑safe multiply into size_t
#if HAVE_UINT128
    __uint128_t prod = ((__uint128_t)base) * ((__uint128_t)mult);
    if (prod > ( __uint128_t)SIZE_MAX) return false;
    *out = (size_t)prod;
#else
    if (base > ULLONG_MAX / mult) return false;
    unsigned long long prod = base * mult;
    if (prod > (unsigned long long)SIZE_MAX) return false;
    *out = (size_t)prod;
#endif
    return true;
}
