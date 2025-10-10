#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dag.h"

/*
 * This program reads a digraph6 code from stdin, one per line, and prints
 * the canonical digraph6 code to stdout.
 */

int main(void) {

    int n;
    char line[1024]; // Buffer for reading input lines
    char d6[128];
    
    if (fgets(line, 1024, stdin) == NULL) {
        return 0;
    }
    n = graphsize(line);
    if (n==0) {
        return 1;
    }
    init_nauty_data(n);
    printf("%s\n", d6_to_d6_canon(line, d6));
    // Loop through each line of standard input
    while (fgets(line, sizeof(line), stdin)) {
        printf("%s\n", d6_to_d6_canon(line, d6)); 
    }

    return 0;
}

