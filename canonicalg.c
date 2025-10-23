#include "common.h"
#include "dag.h"

/*
 * This program reads a digraph6 code from stdin, one per line, and prints
 * the canonical digraph6 code to stdout.
 */

int main(void) {

    int n;
    char line[MAXLINE]; // Buffer for reading input lines
    char d6[MAXLINE];
    
    if (fgets(line, sizeof(line), stdin) == NULL) {
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

