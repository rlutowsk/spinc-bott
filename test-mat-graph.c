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
    // size_t match = 0, unmatch = 0;

    if (fgets(line, sizeof(line), stdin) == NULL) {
        fprintf(stderr, "No input provided\n");
        return 1;
    }

    n = graphsize(line);
    m = SETWORDSNEEDED(n);
    init_nauty_data(n);
    atexit(free_nauty_data);

    graph *g = get_dag_g();
    do {
        // Remove trailing newline
        line[strcspn(line, "\r\n")] = 0;

        // Convert digraph6 string to adjacency matrix
        if (matrix_from_d6(line, mat_from_d6, n)<0) {
            fprintf(stderr, "matrix_from_d6 failed for: %s\n", line);
            exit(1);
        }

        // Encode graph by line
        stringtograph( line, g, m );

        if (matrix_from_graph(g, mat_from_graph) != 0) {
            fprintf(stderr, "matrix_from_graph failed for: %s\n", line);
            exit(1);
        }

        // Compare both keys
        if (memcmp(mat_from_d6, mat_from_graph, n*sizeof(vec_t)) != 0) {
            fprintf(stderr, "conversion mismatch\n");
            exit(1);
        }
        // Convert adjacency matrix back to nauty graph
        matrix_to_graph(g, mat_from_graph);
        // Convert nauty graph back to digraph6 string
        graph_to_d6(g, m, n, buffer);
        if (strcmp(line, buffer) != 0) {
            fprintf(stderr, "stringtograph/graph_to_d6 mismatch\n");
            exit(1);
        }
    } while (fgets(line, sizeof(line), stdin));

    return 0;
}
