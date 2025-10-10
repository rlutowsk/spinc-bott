#include "dag.h"
#include <stdbool.h>
/*
 * DYNALLSTAT in nauty.h is used to declare static variables
 */
DYNALLSTAT(graph, dag_g     , dag_g_sz     );
DYNALLSTAT(graph, dag_canong, dag_canong_sz);
DYNALLSTAT(graph, dag_upperg, dag_upperg_sz);
DYNALLSTAT(int  , dag_lab   , dag_lab_sz   );
DYNALLSTAT(int  , dag_ptn   , dag_ptn_sz   );
DYNALLSTAT(int  , dag_orbits, dag_orbits_sz);
DYNALLSTAT(int  , dag_in_deg, dag_in_deg_sz);
DYNALLSTAT(int  , dag_queue , dag_queue_sz );
DYNALLSTAT(int  , dag_order , dag_order_sz );
DYNALLSTAT(int  , dag_pos   , dag_pos_sz   );

static TLS_ATTR int dag_n;
static TLS_ATTR int dag_m;
static TLS_ATTR int dag_int_arr_sz;

void init_nauty_data(int n)
{
    dag_n = n;
    dag_m = SETWORDSNEEDED(n);

    DYNALLOC1(int  , dag_lab   , dag_lab_sz   , dag_n       , "init_nauty");
    DYNALLOC1(int  , dag_ptn   , dag_ptn_sz   , dag_n       , "init_nauty");
    DYNALLOC1(int  , dag_orbits, dag_orbits_sz, dag_n       , "init_nauty");
    DYNALLOC1(int  , dag_in_deg, dag_in_deg_sz, dag_n       , "init_nauty");
    DYNALLOC1(int  , dag_queue , dag_queue_sz , dag_n       , "init_nauty");
    DYNALLOC1(int  , dag_order , dag_order_sz , dag_n       , "init_nauty");
    DYNALLOC1(int  , dag_pos   , dag_pos_sz   , dag_n       , "init_nauty");
    DYNALLOC2(graph, dag_canong, dag_canong_sz, dag_n, dag_m, "init_nauty");
    DYNALLOC2(graph, dag_g     , dag_g_sz     , dag_n, dag_m, "init_nauty");
    DYNALLOC2(graph, dag_upperg, dag_upperg_sz, dag_n, dag_m, "init_nauty");

    dag_int_arr_sz = dag_n * sizeof(int);
}

void free_nauty_data(void)
{
    DYNFREE(dag_g     , dag_g_sz     );
    DYNFREE(dag_canong, dag_canong_sz);
    DYNFREE(dag_upperg, dag_upperg_sz);
    DYNFREE(dag_lab   , dag_lab_sz   );
    DYNFREE(dag_ptn   , dag_ptn_sz   );
    DYNFREE(dag_orbits, dag_orbits_sz);
    DYNFREE(dag_in_deg, dag_in_deg_sz);
    DYNFREE(dag_queue , dag_queue_sz );
    DYNFREE(dag_order , dag_order_sz );

    dag_int_arr_sz = 0;
}

char* matrix_to_d6(const vec_t *mat, int dim, char *dag_gcode)
{
    int i,j,k;
    char *p,x;
    int n = dim;

    p = dag_gcode;
    *p++ = '&';
    encodegraphsize(n,&p);

    k = 6;
    x = 0;

    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            x <<= 1;
            if (C(mat[j], n, i)) x |= 1;
            if (--k == 0)
            {
                *p++ = (char)(BIAS6 + x);
                k = 6;
                x = 0;
            }
        }
    }

    if (k != 6) *p++ = (char)(BIAS6 + (x << k));

    *p = '\0';

    return dag_gcode;
}

int matrix_from_d6(char *s, vec_t *mat, ind_t dim)
{
    char *p;
    int n,i,j,k,x;

    if (s[0] != '&') {
        return -1;
    }

    n = graphsize(s);
    
    if (n <= 0 || n > dim) {
        return -1;
    }

    memset(mat, 0, dim * sizeof(vec_t)); // Initialize the matrix with zeros

    x = 0;
    p = s + 1 + SIZELEN(n);
    // p = s + SIZELEN(n);

    k = 1;
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            if (--k == 0)
            {
                k = 6;
                x = *(p++) - BIAS6;
            }

            if ((x & TOPBIT6))
            {
                mat[i] |= (1ULL << j);
            }
            x <<= 1;
        }
    }
    return 0;
}

char* graph_to_d6(graph *g, int m, int n, char *dag_gcode)
/* convert nauty graph to digraph6 string, including \0 */
{
    int i,j,k;
    char *p,x;
    set *gj;
    // size_t ii;

    // ii = D6LEN(n)+3;

    // DYNALLOC1(char,dag_gcode,dag_gcode_sz,ii,"graph_to_d6");

    p = dag_gcode;
    *p++ = '&';
    encodegraphsize(n,&p);

    k = 6;
    x = 0;

    for (j = 0; j < n; ++j)
    {
        gj = GRAPHROW(g,j,m);
        for (i = 0; i < n; ++i)
        {
            x <<= 1;
            if (ISELEMENT(gj,i)) x |= 1;
            if (--k == 0)
            {
                *p++ = (char)(BIAS6 + x);
                k = 6;
                x = 0;
            }
        }
    }

    if (k != 6) *p++ = (char)(BIAS6 + (x << k));

    // *p++ = '\n';
    *p = '\0';

    return dag_gcode;
}


/**
 * @brief Generates the canonical labeling of a directed graph.
 *
 * This function sets up the necessary options and memory allocations to compute
 * the canonical form of a directed graph using the nauty library. It prepares
 * label, partition, and orbit arrays, as well as storage for the canonical graph,
 * and then calls `densenauty` to perform the canonical labeling.
 *
 * @param m The number of words per adjacency matrix row (as required by nauty).
 * @param n The number of vertices in the graph.
 *
 * The function uses global variables `dag_*`.
 */
#include <nautinv.h>
static inline void generate_canon_digraph(int m, int n)
{
    DEFAULTOPTIONS_DIGRAPH(dag_options);
    dag_options.getcanon = TRUE;

    statsblk dag_stats;

    densenauty(dag_g, dag_lab, dag_ptn, dag_orbits, &dag_options, &dag_stats, m, n, dag_canong);
}

char* matrix_to_d6_canon(const vec_t *mat, int n, char *dag_gcode)
{
    //int m = SETWORDSNEEDED(n);

    EMPTYGRAPH(dag_g, dag_m, dag_n);

    for (int i=0; i<dag_n; i++) {
        if (row_sum(mat[i]) == 0) {
            continue;
        }
        for (int j=0; j<dag_n; j++) {
            if ( C(mat[i], dag_n, j) ) {
                ADDONEARC(dag_g, i, j, dag_m);
            }
        }
    }
    generate_canon_digraph(dag_m, dag_n);
    return graph_to_d6(dag_canong, dag_m, dag_n, dag_gcode);
}

vec_t* matrix_to_matrix_canon(const vec_t *mat, int n, vec_t *out)
{
    if (out == NULL) {
        return NULL;
    }

    EMPTYGRAPH(dag_g, dag_m, dag_n);

    for (int i=0; i<dag_n; i++) {
        if (row_sum(mat[i]) == 0) {
            continue;
        }
        for (int j=0; j<dag_n; j++) {
            if ( C(mat[i], dag_n, j) ) {
                ADDONEARC(dag_g, i, j, dag_m);
            }
        }
    }
    generate_canon_digraph(dag_m, dag_n);
    
    for (int i=0; i<dag_n; ++i) {
        out[i]  = 0;
        set *gi = GRAPHROW(dag_canong,i,dag_m);
        for (int j=0; j<dag_n; ++j) {
            if (ISELEMENT(gi,j)) {
                out[i] |= 1UL<<j;
            }
        }
    }
    return out;
}

char* d6_to_d6_canon(char *src, char *dst)
{
    EMPTYGRAPH( dag_g, dag_m, dag_n );
    stringtograph( src, dag_g, dag_m );
    generate_canon_digraph( dag_m, dag_n );
    return graph_to_d6( dag_canong, dag_m, dag_n, dst );
}

static bool topo_sort(const graph *g, int n, int m, int *in_degree, int *q, int *order)
{
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
    return outc == n;
}

char *d6_to_d6_upper(char *src, char *dst)
{
    stringtograph( src, dag_g, dag_m );

    memset( dag_in_deg, 0, dag_int_arr_sz );
    memset( dag_queue , 0, dag_int_arr_sz );
    memset( dag_order , 0, dag_int_arr_sz );
    memset( dag_pos   , 0, dag_pos_sz     );

    if (!topo_sort(dag_g, dag_n, dag_m, dag_in_deg, dag_queue, dag_order)) {
        return NULL;
    }
    EMPTYGRAPH( dag_upperg, dag_m, dag_n );
    
    for (int i = 0; i < dag_n; ++i) dag_pos[dag_order[i]] = i;
    
    for (int u = 0; u < dag_n; ++u) {
        set *row = GRAPHROW(dag_g, u, dag_m);
        int iu = dag_pos[u];
        for (int v = 0; v < dag_n; ++v) {
            if (ISELEMENT(row, v)) {
                int iv = dag_pos[v];
                ADDONEARC(dag_upperg, iu, iv, dag_m);
            }
        }
    }
    graph_to_d6(dag_upperg, dag_m, dag_n, dst);
    return dst;
}
