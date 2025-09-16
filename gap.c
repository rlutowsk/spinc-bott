#include "src/compiled.h"
#include "bott.h"

/* default values of dimension, num of threads */
ind_t dim;
vec_t *mat, *aux, *cache;
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
    char buffer[1024];

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
    aux = init(dim);

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

Obj OBJ_MAT(vec_t *mat, ind_t dim)
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

void MAT_OBJ(vec_t *mat, Obj M)
{
    ind_t i, j;
    if (initialized() && IS_PLIST(M) && LEN_PLIST(M)==dim) {
        for (i=0; i<dim; i++) {
            vec_t r = 0;
            Obj row = ELM_PLIST(M, i+1);
            if (IS_PLIST(row) && LEN_PLIST(row)==dim) {
                for (j=0; j<dim; j++) {
                    Obj e = ELM_PLIST(row, j+1);
                    if (IS_INTOBJ(e)) {
                        r |= (INT_INTOBJ(e)&1) << j;
                    } else {
                        ErrorQuit("matrix element is not a small integer", 0L, 0L);
                    }
                }
                mat[i] = r;
            } else {
                ErrorQuit("matrix row is not a list of length %d", (Int)dim, 0L);
            }
        }
    } else {
        ErrorQuit("matrix is not a list of length %d", (Int)dim, 0L);
    }
}

Obj BottMatObj(Obj self)
{
    if (initialized()) {
        return OBJ_MAT(mat, dim);
    }
    return (Obj)0;
}

Obj BottAuxObj(Obj self)
{
    if (initialized()) {
        return OBJ_MAT(aux, dim);
    }
    return (Obj)0;
}

Obj BottSetMat(Obj self, Obj M)
{
    if (initialized()) {
        MAT_OBJ(mat, M);
        return (Obj)0;
    }
    return Fail;
}

Obj BottSetAux(Obj self, Obj M)
{
    if (initialized()) {
        MAT_OBJ(aux, M);
        return (Obj)0;
    }
    return Fail;
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
    mat[dim-2] = 0;
    mat[dim-1] = 0;
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

Obj BottLastNonzeroRow(Obj self)
{
    ind_t c;
    for (c=dim-1; c>=0; c--) {
        if (mat[c]>0) {
            return INTOBJ_INT(c+1);
        }
    }
    return Fail;
}

Obj BottSwapRowsAndCols(Obj self, Obj r1, Obj r2)
{
    ind_t row1 = INT_INTOBJ(r1);
    ind_t row2 = INT_INTOBJ(r2);
    if (uninitialized() || row1<1 || row1>dim || row2<1 || row2>dim) {
        return Fail;
    }
    swap_rows_and_cols(mat, aux, dim, row1-1, row2-1);
    return (Obj)0;
}

Obj BottConditionalAddCol(Obj self, Obj k)
{
    ind_t col = INT_INTOBJ(k);
    if (uninitialized() || col<1 || col>dim) {
        return Fail;
    }
    conditional_add_col(mat, aux, dim, col-1);
    return (Obj)0;
}

Obj BottConditionalAddRow(Obj self, Obj l, Obj m)
{
    ind_t row_l = INT_INTOBJ(l);
    ind_t row_m = INT_INTOBJ(m);
    if (uninitialized() || row_l<1 || row_l>dim || row_m<1 || row_m>dim) {
        return Fail;
    }
    conditional_add_row(mat, aux, dim, row_l-1, row_m-1);
    return (Obj)0;
}

Obj BottEncodeMat(Obj self)
{
    static char buffer[1024];
    if (initialized()) {
        encode_matrix(mat, dim, buffer);
        return MakeString(buffer);
    }
    return Fail;
}

Obj BottEncodeAux(Obj self)
{
    static char buffer[1024];
    if (initialized()) {
        encode_matrix(aux, dim, buffer);
        return MakeString(buffer);
    }
    return Fail;
}

Obj BottDecodeMat(Obj self, Obj code)
{
    if (initialized()) {
        decode_matrix(mat, dim, CSTR_STRING(code));
    }
    return (Obj)0;
}

int is_upper_triangular(const vec_t *mat, const ind_t dim)
{
    vec_t row_mask = 1;
    for (ind_t i=0; i<dim; i++) {
        row_mask |= (1<<i);
        if ( mat[i] & row_mask ) {
            return 0;
        }
    }
    return 1;
}

Obj BottIsUpperTriangularMat(Obj self)
{
    return is_upper_triangular(mat, dim) ? True : False;
}

Obj BottIsUpperTriangularAux(Obj self)
{
    return is_upper_triangular(aux, dim) ? True : False;
}

Obj BottIsUpperTriangularByCode(Obj self, Obj code)
{
    decode_matrix(aux, dim, CSTR_STRING(code));
    return is_upper_triangular(aux, dim) ? True : False;
}

/* 
 * GVarFunc - list of functions to export
 */
static StructGVarFunc GVarFunc[] = {
    { "BottInit"          , 1, "dim"   , BottInit           , "bott.c:BottInit" },
    { "BottClear"         , 0, ""      , BottClear          , "bott.c:BottClear" },
    { "BottPrint"         , 0, ""      , BottPrint          , "bott.c:BottPrint" },
    { "BottGetMat"        , 0, ""      , BottMatObj         , "bott.c:BottGetMat" },
    { "BottSetMat"        , 1, "M"     , BottSetMat         , "bott.c:BottSetMat" },
    { "BottGetAux"        , 0, ""      , BottAuxObj         , "bott.c:BottGetAux" },
    { "BottSetAux"        , 1, "M"     , BottSetAux         , "bott.c:BottSetAux" },
    { "BottNext"          , 0, ""      , BottNext           , "bott.c:BottNext" },
    { "BottSetState"      , 1, "state" , BottSetState       , "bott.c:BottSetState" },
    { "BottGetState"      , 0, ""      , BottGetState       , "bott.c:BottGetState" },
    { "BottIsSpinc"       , 0, ""      , BottIsSpinc        , "bott.c:BottIsSpinc" },
    { "BottIsSpin"        , 0, ""      , BottIsSpin         , "bott.c:BottIsSpin" },
    { "BottLastNonzeroRow", 0, ""      , BottLastNonzeroRow , "bott.c:BottLastNonzeroRow" },
    { "BottSwapRowsAndCols", 2, "r1,r2", BottSwapRowsAndCols, "bott.c:BottSwapRowsAndCols" },
    { "BottConditionalAddCol", 1, "k", BottConditionalAddCol, "bott.c:BottConditionalAddCol" },
    { "BottConditionalAddRow", 2, "l,m", BottConditionalAddRow, "bott.c:BottConditionalAddRow" },
    { "BottEncodeMat"     , 0, ""      , BottEncodeMat      , "bott.c:BottEncodeMat" },
    { "BottEncodeAux"     , 0, ""      , BottEncodeAux      , "bott.c:BottEncodeAux" },
    { "BottDecodeMat"     , 1, "code"  , BottDecodeMat      , "bott.c:BottDecodeMat" },
    { "BottIsUpperTriangularMat", 0, "", BottIsUpperTriangularMat, "bott.c:BottIsUpperTriangularMat" },
    { "BottIsUpperTriangularAux", 0, "", BottIsUpperTriangularAux, "bott.c:BottIsUpperTriangularAux" },
    { "BottIsUpperTriangularByCode", 1, "code", BottIsUpperTriangularByCode, "bott.c:BottIsUpperTriangularByCode" },
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
