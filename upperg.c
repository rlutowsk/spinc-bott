#include "dag.h"

int main(void)
{
    int   n;
    char  s[MAXLINE];
    char d6[MAXLINE];

    if ( fgets(s, sizeof(s), stdin)==NULL ) {
        fprintf(stderr, "file read error, quitting...\n");
        exit(1);
    }
    n = graphsize(s);
    if (n <= 0) {
        fprintf(stderr, "read non-positive number of vertices, quitting...\n");
        exit(1);
    }
    init_nauty_data(n);

    do {
        if (d6_to_d6_upper(s, d6)==NULL) {
            fprintf(stderr, "%s: error in converting graph to topological order, quitting...\n", s);
            free_nauty_data();
            exit(1);
        }
        printf("%s\n", d6);
    } while (fgets(s, sizeof(s), stdin) != NULL);

    free_nauty_data();

    return 0;
}
