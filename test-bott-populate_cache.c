#include "common.h"
#include "bott.h"
#include "dag.h"
#include "bucket.h"
#include "adjpack11.h"

int main(void)
{
    const int    dim =  8;
    const size_t siz = 64;
    const vec_t  vectors[] = {
        /* dimension 1 */
        0x00, //00000000
        0xc0, //11000000
        /* dimension 2 */
        0x60, //01100000
        0xa0, //10100000
        /* dimension 3 */
        0x30, //00110000
        0xf0, //11110000
        0x50, //01010000
        0x90, //10010000
        /* dimension 4 */
        0x18, //00011000
        0xd8, //11011000
        0x78, //01111000
        0xb8, //10111000
        0x28, //00101000
        0xe8, //11101000
        0x48, //01001000
        0x88, //10001000
        /* dimension 5 */
        0x0c, //00001100
        0xcc, //11001100
        0x6c, //01101100
        0xac, //10101100
        0x3c, //00111100
        0xfc, //11111100
        0x5c, //01011100
        0x9c, //10011100
        0x14, //00010100
        0xd4, //11010100
        0x74, //01110100
        0xb4, //10110100
        0x24, //00100100
        0xe4, //11100100
        0x44, //01000100
        0x84, //10000100
        /* dimensions 6 */
        0x06, //00000110
        0xc6, //11000110
        0x66, //01100110
        0xa6, //10100110
        0x36, //00110110
        0xf6, //11110110
        0x56, //01010110
        0x96, //10010110
        0x1e, //00011110
        0xde, //11011110
        0x7e, //01111110
        0xbe, //10111110
        0x2e, //00101110
        0xee, //11101110
        0x4e, //01001110
        0x8e, //10001110
        0x0a, //00001010
        0xca, //11001010
        0x6a, //01101010
        0xaa, //10101010
        0x3a, //00111010
        0xfa, //11111010
        0x5a, //01011010
        0x9a, //10011010
        0x12, //00010010
        0xd2, //11010010
        0x72, //01110010
        0xb2, //10110010
        0x22, //00100010
        0xe2, //11100010
        0x42, //01000010
        0x82  //10000010 
    };
    

    vec_t *cache = NULL;
    size_t cache_size = 0;
    printf("=== [bott-populate_cache] testing dimension %d ===\n", dim);
    populate_cache(&cache, &cache_size, dim);
    if (cache_size != siz) {
        fprintf(stderr, "    Size mismatch. Expected: %lu, calculated: %lu.\n", siz, cache_size);
        exit(1);
    }
    for (size_t i=0; i<cache_size; ++i) {
        if (cache[i] != vectors[i]) {
            fprintf(stderr, "    Cache miscalculations at %lu. Expected: %lx, calculated: %lx.\n", i, vectors[i], cache[i]);
            exit(1);
        }
    }
    free(cache);

    printf("=== [bott-populate_cache] all tests passed ===\n");
    
    return 0;
}
