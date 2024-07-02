#include "src/compiled.h"
#include "bott.h"

/* default values of dimension, num of threads */
ind_t dim;
vec_t *mat, *cache;
/* state is a number which is used to generate RBM matrix */
vec_t state, max_state;
size_t spinc, spin;
size_t cache_size;

inline int uninitialized(void)
{
    return (cache_size==0);
}

inline int initialized(void)
{
    return (cache_size>0);
}

Obj BottInit(Obj self, Obj d)
{
    ind_t c;

    if (initialized()) {
        return False;
    }
    if (! IS_INTOBJ(d)) {
        ErrorQuit("dim must be a small integer (not a %s)", (Int)TNAM_OBJ(d), 0L);
    }
    dim = INT_INTOBJ(d);
    if (dim<3 || dim>10) {
        ErrorQuit("value of dim must be between 3 and 10", 0L, 0L);
    }

    for (c=0, max_state=1; c<dim-1; c++) {
        max_state <<= (dim-c-2);
    }
    max_state -= 1;

    populate_cache(&cache, &cache_size, dim);

    mat = init(dim);

    state = 0;

    return True;
}

Obj BottClear(Obj self)
{
    if (initialized()) {
        free(cache);
        free(mat);
        cache_size = 0;
    }
    return (Obj)0;
}

Obj BottPrint(Obj self)
{
    ind_t i, j;
    const char c[] = { '.', '1' };

    if (initialized()) {
        for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {
                printf("%c", c[C(mat[i],dim,j)]);
            }
            printf("\n");
        }
    }
    return (Obj)0;
}

Obj BottMat(Obj self)
{
    ind_t i,j;
    Obj M, row;

    M = NEW_PLIST(T_PLIST, dim);
    SET_LEN_PLIST(M, dim);
    row = NEW_PLIST(T_PLIST, dim);
    SET_LEN_PLIST(row, dim);
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            SET_ELM_PLIST( row, j+1, INTOBJ_INT(C(mat[i],dim,j)));
        }
        SET_ELM_PLIST(M, i+1, SHALLOW_COPY_OBJ(row) );
        CHANGED_BAG(M);
    }
    return M;
}

Obj BottNext(Obj self)
{
    if (uninitialized() || state==max_state) {
        return Fail;
    }
    set(mat, cache, ++state, dim);
    return INTOBJ_INT(state);
}

Obj BottSetState(Obj self, Obj s)
{
    vec_t candidate = INT_INTOBJ(s);
    if (uninitialized() || candidate<0 || candidate>=max_state) {
        return Fail;
    }
    state = candidate;
    set(mat, cache, state, dim);
    return s;
}

Obj BottGetState(Obj self)
{
    if (uninitialized()) {
        return Fail;
    }
    return INTOBJ_INT(state);
}

Obj BottIsSpinc(Obj self)
{
    return (is_spinc(mat, dim))?True:False;
}

Obj BottIsSpin(Obj self)
{
    return (is_spin(mat, dim))?True:False;
}

/* 
 * GVarFunc - list of functions to export
 */
static StructGVarFunc GVarFunc[] = {
    { "BottInit" ,     1, "dim"   , BottInit     , "bott.c:BottInit"     },
    { "BottClear",     0, ""          , BottClear    , "bott.c:BottClear"    },
    { "BottPrint",     0, ""          , BottPrint    , "bott.c:BottPrint"    },
    { "BottMat"  ,     0, ""          , BottMat      , "bott.c:BottMat"      },
    { "BottNext" ,     0, ""          , BottNext     , "bott.c:BottNext"     },
    { "BottSetState" , 1, "state"     , BottSetState , "bott.c:BottSetState" },
    { "BottGetState" , 0, ""          , BottGetState , "bott.c:BottGetState" },
    { "BottIsSpinc"  , 0, ""          , BottIsSpinc  , "bott.c:BottIsSpinc"  },
    { "BottIsSpin"   , 0, ""          , BottIsSpin   , "bott.c:BottIsSpin"   },
    { 0 }
};

static Int InitKernel (StructInitInfo * module)
{
    InitHdlrFuncsFromTable(GVarFunc);
    return 0;
}

static Int InitLibrary(StructInitInfo * module)
{
    InitGVarFuncsFromTable(GVarFunc);

    mat   = NULL;
    cache = NULL;
    dim   = 6;
    cache_size = 0;

    return 0;
}

static StructInitInfo module = {
    MODULE_DYNAMIC,
    "shared",
    0,
    0,
    0,
    0,
    InitKernel,
    InitLibrary,
    0,
    0,
    0,
    0
};

StructInitInfo * Init__Dynamic(void)
{
  return &module;
}
