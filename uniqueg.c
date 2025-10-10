#include <stdio.h>
#include <stdbool.h>
#include <glib.h>
#include "dag.h"


int main(void) {
    char s[1024];

    // create hash table
    GHashTable* code_set = g_hash_table_new_full(
        g_str_hash,    // Hash function for keys (strings)
        g_str_equal,   // Equality function for keys (strings)
        g_free,        // Function to free the key when the table is destroyed
        NULL           // Function to free the value (not needed for simple booleans)
    );

    // Main loop to read one line (one graph) at a time
    while (fgets(s, sizeof(s), stdin) != NULL) {
        if ( g_hash_table_insert(code_set, g_strdup(s), GINT_TO_POINTER(TRUE))==TRUE ) {
            printf("%s", s);
        }
    }
    // clean hash table
    g_hash_table_unref(code_set); 
    return 0;
}
