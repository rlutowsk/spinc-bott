#include "common.h"
#include "dag.h"

// Topologiczne sortowanie (Kahn), działa dla grafów skierowanych zapisanych
// w g (n wierzch., m = SETWORDSNEEDED(n)). Zwraca true dla DAG.
// Wynikowa permutacja wierzchołków w tablicy 'order' (order[i] = oryginalny wierzchołek na miejscu i)
/*
static bool topo_sort(const graph *g, int n, int m, int *order)
{
    int *in_degree = (int*)calloc((size_t)n, sizeof(int));
    int *q = (int*)malloc((size_t)n * sizeof(int));
    if (!in_degree || !q) { free(in_degree); free(q); return false; }

    // 1) stopnie wejściowe
    for (int u = 0; u < n; ++u) {
        set *row = GRAPHROW(g, u, m);
        for (int v = 0; v < n; ++v) {
            if (ISELEMENT(row, v)) {
                ++in_degree[v];
            }
        }
    }

    // 2) kolejka wierzchołków o stopniu 0
    int head = 0, tail = 0;
    for (int v = 0; v < n; ++v) {
        if (in_degree[v] == 0) q[tail++] = v;
    }

    // 3) przetwarzanie
    int outc = 0;
    while (head < tail) {
        int u = q[head++];
        order[outc++] = u;

        set *row = GRAPHROW(g, u, m);
        for (int v = 0; v < n; ++v) {
            if (ISELEMENT(row, v)) {
                if (--in_degree[v] == 0) q[tail++] = v;
            }
        }
    }

    free(in_degree);
    free(q);
    return outc == n;
}
*/
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
