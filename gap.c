#include "src/compiled.h"
#include "bott.h"
#include "dag.h"

/* default values of dimension, num of threads */
static ind_t dim;
static vec_t *mat, *aux, *cache;
/* state is a number which is used to generate RBM matrix */
static vec_t state, max_state;
static size_t spinc, spin;
static size_t cache_size;

static char d6[128];

Obj code_mat, code_aux;
size_t code_len;

static inline int uninitialized(void)
{
    return (cache_size==0);
}

static inline int initialized(void)
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

    code_len = mat_characters[dim]; // +1 for '\0'

    code_mat = NEW_STRING(code_len);
    CHARS_STRING(code_mat)[code_len] = '\0';

    code_aux = NEW_STRING(code_len);
    CHARS_STRING(code_aux)[code_len] = '\0';

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

Obj PLIST_MAT(vec_t *mat, ind_t dim)
{
    ind_t i,j,k;
    int len;
    Obj M, row;

    M = NEW_PLIST(T_PLIST, dim);
    SET_LEN_PLIST(M, dim);

    for (i=0; i<dim; ) {
        len = row_sum(mat[i]);
        if (len>0) {
            row = NEW_PLIST(T_PLIST, len);
            SET_LEN_PLIST(row, len);
        } else {
            row = NewEmptyPlist();
        }

        for (j=0, k=0; k<len && j<dim; ) {
            // increment j after checking
            if (C(mat[i],dim,j++)) {
                // if here, the j value is already incremented
                SET_ELM_PLIST( row, ++k, INTOBJ_INT(j));
            }
        }
        // i is incremented here
        SET_ELM_PLIST(M, ++i, row );
        CHANGED_BAG(M);
    }
    return M;
}

Obj BottListMat(Obj self)
{
    if (initialized()) {
        return PLIST_MAT(mat, dim);
    }
    return (Obj)0;
}

Obj BottListAux(Obj self)
{
    if (initialized()) {
        return PLIST_MAT(aux, dim);
    }
    return (Obj)0;
}

Obj OBJ_MAT(vec_t *mat, ind_t dim)
{
    ind_t i,j;
    Obj M, row;

    M = NEW_PLIST(T_PLIST, dim);
    SET_LEN_PLIST(M, dim);

    for (i=0; i<dim; i++) {
        row = NEW_PLIST(T_PLIST, dim);
        SET_LEN_PLIST(row, dim);

        for (j=0; j<dim; j++) {
            SET_ELM_PLIST( row, j+1, INTOBJ_INT(C(mat[i],dim,j)));
        }
        SET_ELM_PLIST(M, i+1, row );
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
    matrix_by_state(mat, cache, ++state, dim);
    return INTOBJ_INT(state);
}

Obj BottSetState(Obj self, Obj s)
{
    vec_t candidate = INT_INTOBJ(s);
    if (uninitialized() || candidate<0 || candidate>=max_state) {
        return Fail;
    }
    state = candidate;
    matrix_by_state(mat, cache, state, dim);
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

static inline Obj are_not_equal()
{
    for (ind_t i=0; i<dim; i++) {
        if ( mat[i] != aux[i] ) {
            return True;
        }
    }
    return False;
}

Obj BottSwapRowsAndCols(Obj self, Obj r1, Obj r2)
{
    ind_t row1 = INT_INTOBJ(r1);
    ind_t row2 = INT_INTOBJ(r2);
    if (uninitialized() || row1<1 || row1>dim || row2<1 || row2>dim) {
        return Fail;
    }
    swap_rows_and_cols(mat, aux, dim, row1-1, row2-1);
    //return are_not_equal(); 
    return (Obj)0;
}

Obj BottSwapRowsAndColsAux(Obj self, Obj r1, Obj r2)
{
    ind_t row1 = INT_INTOBJ(r1);
    ind_t row2 = INT_INTOBJ(r2);
    if (uninitialized() || row1<1 || row1>dim || row2<1 || row2>dim) {
        return Fail;
    }
    swap_rows_and_cols(aux, aux, dim, row1-1, row2-1);
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
    if ( conditional_add_row(mat, aux, dim, row_l-1, row_m-1) ) {
        return True;
    }
    return False;
}

Obj BottEncodeMat(Obj self)
{
    char *dst;
    if (initialized()) {
        dst = CHARS_STRING(code_mat);
        encode_matrix(mat, dim, dst);
        return code_mat;
    }
    return Fail;
}

Obj BottEncodeAux(Obj self)
{
    static char buffer[1024];
    char *dst;
    if (initialized()) {
        dst = CHARS_STRING(code_aux);
        encode_matrix(aux, dim, dst);
        return code_aux;
    }
    return Fail;
}

Obj BottCopyMatCode(Obj self)
{
    return SHALLOW_COPY_OBJ(code_mat);
    //return MakeStringWithLen( CHARS_STRING(code_mat), code_len );
}

Obj BottCopyAuxCode(Obj self)
{
    return SHALLOW_COPY_OBJ(code_aux);
    //return MakeStringWithLen( CHARS_STRING(code_aux), code_len );
}

Obj BottDecodeMat(Obj self, Obj code)
{
    if (initialized()) {
        decode_matrix(mat, dim, CSTR_STRING(code));
    }
    return (Obj)0;
}

Obj BottDigraph6Mat(Obj self)
{
    if (initialized()) {
        return MakeString( matrix_to_digraph6(mat, dim) );
    }
    return Fail;
}

Obj BottDigraph6Aux(Obj self)
{
    if (initialized()) {
        return MakeString( matrix_to_digraph6(aux, dim) );
    }
    return Fail;
}

Obj BottCanonicalDigraph6Mat(Obj self)
{
    if (initialized()) {
        return MakeString( matrix_to_d6_canon(mat, dim, d6) );
    }
    return Fail;
}

Obj BottCanonicalDigraph6Aux(Obj self)
{
    if (initialized()) {
        return MakeString( matrix_to_d6_canon(aux, dim, d6) );
    }
    return Fail;
}

Obj BottSetMatByDigraph6(Obj self, Obj s)
{
    if (initialized() && IS_STRING(s) && GET_LEN_STRING(s)>2) {
        if ( digraph6_to_matrix(CSTR_STRING(s), mat, dim) == 0 ) {
            return True;
        }
    }
    return Fail;
}

Obj BottSetAuxByDigraph6(Obj self, Obj s)
{
    if (initialized() && IS_STRING(s) && GET_LEN_STRING(s)>2) {
        if ( digraph6_to_matrix(CSTR_STRING(s), aux, dim) == 0 ) {
            return True;
        }
    }
    return Fail;
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

Obj BottWeightMat(Obj self)
{
    if (initialized()) {
        return INTOBJ_INT(matrix_weight(mat, dim));
    }
    return Fail;
}

Obj BottWeightAux(Obj self)
{
    if (initialized()) {
        return INTOBJ_INT(matrix_weight(aux, dim));
    }
    return Fail;
}

Obj BottDimByDigraph6(Obj self, Obj str)
{
    return INTOBJ_INT( graphsize(CSTR_STRING(str)) );
}

/* 
 * GVarFunc - list of functions to export
 */
static StructGVarFunc GVarFunc[] = {
    { "BottInit"                   , 1, "dim"   , BottInit                   , "bott.c:BottInit" },
    { "BottClear"                  , 0, ""      , BottClear                  , "bott.c:BottClear" },
    { "BottPrint"                  , 0, ""      , BottPrint                  , "bott.c:BottPrint" },
    { "BottListMat"                , 0, ""      , BottListMat                , "bott.c:BottListMat" },
    { "BottListAux"                , 0, ""      , BottListAux                , "bott.c:BottListAux" },
    { "BottGetMat"                 , 0, ""      , BottMatObj                 , "bott.c:BottGetMat" },
    { "BottSetMat"                 , 1, "M"     , BottSetMat                 , "bott.c:BottSetMat" },
    { "BottGetAux"                 , 0, ""      , BottAuxObj                 , "bott.c:BottGetAux" },
    { "BottSetAux"                 , 1, "M"     , BottSetAux                 , "bott.c:BottSetAux" },
    { "BottNext"                   , 0, ""      , BottNext                   , "bott.c:BottNext" },
    { "BottSetState"               , 1, "state" , BottSetState               , "bott.c:BottSetState" },
    { "BottGetState"               , 0, ""      , BottGetState               , "bott.c:BottGetState" },
    { "BottIsSpinc"                , 0, ""      , BottIsSpinc                , "bott.c:BottIsSpinc" },
    { "BottIsSpin"                 , 0, ""      , BottIsSpin                 , "bott.c:BottIsSpin" },
    { "BottLastNonzeroRow"         , 0, ""      , BottLastNonzeroRow         , "bott.c:BottLastNonzeroRow" },
    { "BottSwapRowsAndCols"        , 2, "r1,r2" , BottSwapRowsAndCols        , "bott.c:BottSwapRowsAndCols" },
    { "BottSwapRowsAndColsAux"     , 2, "r1,r2" , BottSwapRowsAndColsAux     , "bott.c:BottSwapRowsAndColsAux" },
    { "BottConditionalAddCol"      , 1, "k"     , BottConditionalAddCol      , "bott.c:BottConditionalAddCol" },
    { "BottConditionalAddRow"      , 2, "l,m"   , BottConditionalAddRow      , "bott.c:BottConditionalAddRow" },
    { "BottEncodeMat"              , 0, ""      , BottEncodeMat              , "bott.c:BottEncodeMat" },
    { "BottEncodeAux"              , 0, ""      , BottEncodeAux              , "bott.c:BottEncodeAux" },
    { "BottDigraph6Mat"            , 0, ""      , BottDigraph6Mat            , "bott.c:BottDigraph6Mat" },
    { "BottDigraph6Aux"            , 0, ""      , BottDigraph6Aux            , "bott.c:BottDigraph6Aux" },
    { "BottCanonicalDigraph6Mat"   , 0, ""      , BottCanonicalDigraph6Mat   , "bott.c:BottCanonicalDigraph6Mat" },
    { "BottCanonicalDigraph6Aux"   , 0, ""      , BottCanonicalDigraph6Aux   , "bott.c:BottCanonicalDigraph6Aux" },
    { "BottSetMatByDigraph6"       , 1, "s"     , BottSetMatByDigraph6       , "bott.c:BottSetMatByDigraph6" },
    { "BottSetAuxByDigraph6"       , 1, "s"     , BottSetAuxByDigraph6       , "bott.c:BottSetAuxByDigraph6" },
    { "BottCopyMatCode"            , 0, ""      , BottCopyMatCode            , "bott.c:BottCopyMatCode" },
    { "BottCopyAuxCode"            , 0, ""      , BottCopyAuxCode            , "bott.c:BottCopyAuxCode" },
    { "BottDecodeMat"              , 1, "code"  , BottDecodeMat              , "bott.c:BottDecodeMat" },
    { "BottIsUpperTriangularMat"   , 0, ""      , BottIsUpperTriangularMat   , "bott.c:BottIsUpperTriangularMat" },
    { "BottIsUpperTriangularAux"   , 0, ""      , BottIsUpperTriangularAux   , "bott.c:BottIsUpperTriangularAux" },
    { "BottIsUpperTriangularByCode", 1, "code"  , BottIsUpperTriangularByCode, "bott.c:BottIsUpperTriangularByCode" },
    { "BottWeightMat"              , 0, ""      , BottWeightMat              , "bott.c:BottWeightMat" },
    { "BottWeightAux"              , 0, ""      , BottWeightAux              , "bott.c:BottWeightAux" },
    { "BottDimByDigraph6"          , 1, "code"  , BottDimByDigraph6          , "bott.c:BottDimByDigraph6" },
    { 0 }
};

static Int InitKernel (StructInitInfo * module)
{
    InitHdlrFuncsFromTable(GVarFunc);

    InitGlobalBag(&code_mat, "rbm/gap.c:BOTT_CODE_MAT");
    InitGlobalBag(&code_aux, "rbm/gap.c:BOTT_CODE_AUX");

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
