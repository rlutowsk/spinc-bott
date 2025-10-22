#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include "dag.h"
#include "adjpack11.h"
#include "d6pack11.h"  // for d6pack_decode

int main(void) {
    char line[128], buffer[128];

    vec_t mat_from_d6[11], mat_from_graph[11]; // max n = 11
    int n, m;

    n = 1;
    m = SETWORDSNEEDED(n);
    init_nauty_data(n);
    atexit(free_nauty_data);

    graph *g = get_dag_g();
    FILE *in = fopen("test.d6", "r");
    if (in == NULL) {
        in = stdin;
        printf("=== [mat-graph] testing digraph6 strings from stdin ===\n");
    } else {
        printf("=== [mat-graph] testing digraph6 strings from file 'test.d6' ===\n");
    }
    fflush(stdout);

    while (fgets(line, sizeof(line), in)) {
        if (n != graphsize(line)) {
            n = graphsize(line);
            m = SETWORDSNEEDED(n);
            free_nauty_data();
            init_nauty_data(n);
            g = get_dag_g();
        }

        // Convert digraph6 string to adjacency matrix
        if (matrix_from_d6(line, mat_from_d6, n)<0) {
            fprintf(stderr, "    matrix_from_d6 failed for: %s\n", line);
            exit(1);
        }

        // Encode graph by line
        stringtograph( line, g, m );

        if (matrix_from_graph(g, mat_from_graph) != 0) {
            fprintf(stderr, "    matrix_from_graph failed for: %s\n", line);
            exit(1);
        }

        // Compare both keys
        if (memcmp(mat_from_d6, mat_from_graph, n*sizeof(vec_t)) != 0) {
            fprintf(stderr, "    conversion mismatch\n");
            exit(1);
        }
        // Convert adjacency matrix back to nauty graph
        matrix_to_graph(g, mat_from_graph);
        // Convert nauty graph back to digraph6 string
        graph_to_d6(g, m, n, buffer);
        if (strncmp(line, buffer, strlen(buffer)) != 0) {
            fprintf(stderr, "    stringtograph/graph_to_d6 mismatch\n");
            fprintf(stderr, "      original:   %s", line);
            fprintf(stderr, "      calculated: %s\n", buffer);
            exit(1);
        }
    }
    printf("=== [mat-graph] all tests passed. ===\n");
    if (in != stdin) {
        fclose(in);
    }
    return 0;
}
