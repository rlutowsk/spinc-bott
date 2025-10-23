#include "dag.h"
#include "adjpack11.h"
#include "d6pack11.h"

_Static_assert(sizeof(vec_t) == sizeof(setword), "vec_t and setword types must match for memcpy optimization");
// check at compile time if everything is as expected
_Static_assert(sizeof(vec_t) == 8, "vec_t must be a 64-bit type.");
_Static_assert(sizeof(set)   == 8, "nauty's set type must be 64-bit.");

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

static TLS_ATTR int   dag_n;
static TLS_ATTR int   dag_m;
static TLS_ATTR int   dag_int_arr_sz;
static TLS_ATTR vec_t dag_mask;

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

    dag_mask = (dag_n == 64) ? (vec_t)-1 : (((vec_t)1 << dag_n) - 1);

    // Assertion checked at runtime (in debug mode)
    // Documents a key limitation of the algorithm.
    assert(dag_m == 1 && "This function is optimized for n <= 64 and will not work correctly otherwise.");
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
    DYNFREE(dag_pos   , dag_pos_sz   );

    dag_int_arr_sz = 0;

    dag_mask = 0;
    dag_n = 0;
    dag_m = 0;
}

graph *get_dag_g(void)
{
    return dag_g;
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
            // if (C(mat[j], n, i)) x |= 1;
            x |= C(mat[j], n, i);
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

int matrix_from_d6(char *s, vec_t * __restrict mat, ind_t dim)
{
    key128_t k;
    unsigned n = 0;
    d6pack_decode(s, &k, &n);
    if (n == 0 || (ind_t)n > dim) return -1;
    adjpack_to_matrix(&k, mat, n);
    return 0;
}

/*
#include <string.h>
#include <stdint.h>

// Założenia jak u Ciebie:
// - graphsize(), SIZELEN(n)
// - BIAS6 (63), TOPBIT6 (0x20)
// - typy: vec_t (np. uint64_t) i ind_t

int matrix_from_d6(char *s, vec_t * __restrict mat, ind_t dim)
{
    if (!s || s[0] != '&') return -1;

    const int n = graphsize(s);
    if (n <= 0 || (ind_t)n > dim) return -1;

    // Zerowanie tylko potrzebnej części (lub całego dim – jak wolisz)
    memset(mat, 0, (size_t)dim * sizeof(vec_t));

    const unsigned char *p = (const unsigned char*)s + 1 + SIZELEN(n);

    // Globalny, ciągły strumień bitów: k = ile bitów pozostało w x (0..6)
    unsigned k = 0;
    unsigned x = 0;

    for (int i = 0; i < n; ++i) {
        vec_t row  = 0;
        vec_t mask = 1;

        for (int j = 0; j < n; ++j) {
            if (k == 0) {                 // doładuj kolejne 6 bitów
                x = (unsigned)(*p++) - (unsigned)BIAS6;
                k = 6;
            }
            if (x & TOPBIT6) row |= mask; // testuj MSB (bit5), jak w Twoim kodzie
            x <<= 1;                       // przejdź do kolejnego bitu: bit4 -> bit5 ...
            --k;                           // jeden bit mniej do wykorzystania
            mask <<= 1;                    // następna kolumna w wierszu
        }
        mat[i] = row;                      // jeden zapis na wiersz
    }

    return 0;
}
*/

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
    /* generate graph from mat */
    matrix_to_graph(dag_g, mat);
    /* canonize the graph */
    generate_canon_digraph(dag_m, dag_n);
    return graph_to_d6(dag_canong, dag_m, dag_n, dag_gcode);
}

vec_t* matrix_to_matrix_canon(const vec_t *mat, int n, vec_t *out)
{
    if (out == NULL) {
        return NULL;
    }
    /* generate graph from mat */
    matrix_to_graph(dag_g, mat);
    /* canonize the graph */
    generate_canon_digraph(dag_m, dag_n);
    /* generate out from canonical form of the graph */
    matrix_from_graph(dag_canong, out);
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
    // 1) input degrees
    for (int u = 0; u < n; ++u) {
        set *row = GRAPHROW(g, u, m);
        for (int v = 0; v < n; ++v) {
            if (ISELEMENT(row, v)) {
                ++in_degree[v];
            }
        }
    }

    // 2) queue of vertices with degree 0
    int head = 0, tail = 0;
    for (int v = 0; v < n; ++v) {
        if (in_degree[v] == 0) q[tail++] = v;
    }

    // 3) processing
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

int matrix_from_graph(graph *g, vec_t *mat)
{
    if (g == NULL || mat == NULL) {
        return -1;
    }

    for (int i = 0; i < dag_n; ++i) {
        set *gi = GRAPHROW(g,i,dag_m);
        mat[i] = bitreverse64((uint64_t)*gi) & dag_mask;
    }
    return 0;
}

int matrix_to_graph(graph *g, const vec_t *mat)
{
    if (g == NULL || mat == NULL) {
        return -1;
    }

    for (int i = 0; i < dag_n; ++i) {
        set *gi = GRAPHROW(g,i,dag_m);
        *gi = bitreverse64((uint64_t)(mat[i] & dag_mask));
    }
    return 0;
}