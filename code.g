BottAVLInit := function( dim, info... )
    local t, s;

    BottClear();
    BottInit( dim );

    if IsBound( info[1] ) then
        SetInfoLevel( BottCalcInfo, info[1] );
    fi;

    s := 2^((dim-1)*(dim-2)/2);

    Info( BottCalcInfo, 1, "Populating tree with size ", s, " ... ");
    t := AVLTree( rec( allocsize := s ) );
    repeat 
        BottEncodeMat(); 
        AVLAdd( t, BottCopyMatCode(), true ); 
    until BottNext()=fail;
    Info( BottCalcInfo, 1, "Population completed.");

    BottSetState( 0 );

    return t;
end;

BottAVLNextDel := function(t, dim)
    local orbit, x, i, j, code, seen, rval, first, tri, pos, weight, w, len, AddCode;

    AddCode := function()
        local s, c;
        s := BottEncodeAux();
        if not IsBound(seen.(s)) then
            c := BottCopyAuxCode();
            Add(orbit, c);
            seen.(c) := true;
            if BottIsUpperTriangularAux() then
                AVLDelete( t, c );
            fi;
        fi;
    end;

    code := AVLIndexDelete( t, 1 );
    if code = fail then
        return fail;
    fi;
    orbit:= [ code ];

    seen := rec();
    seen.(code) := true;

    while orbit<>[] do
        BottDecodeMat(Remove(orbit));

        i := 0;
        while i<dim do
            i := i+1;
            j := i;
            while j<dim do
                j := j+1;
                BottSwapRowsAndCols(i, j);
                AddCode();
            od;
        od;

        i := 0;
        while i<dim do
            i := i+1;
            BottConditionalAddCol(i);
            AddCode();
        od;

        i := 0;
        while i<dim do
            i := i+1;
            j := i;
            while j < dim do
                j := j+1;
                if BottConditionalAddRow(i, j) then
                    AddCode();
                fi;
                if BottConditionalAddRow(j, i) then
                    AddCode();
                fi;
            od;
        od;
    od;

    return code;
end;

BottAVLNextPos := function(t, dim)
    local orbit, i, j, code, seen, pos, AddCode;

    AddCode := function()
        local s, c;
        s := BottEncodeAux();
        if not IsBound(seen.(s)) then
            c := BottCopyAuxCode();
            Add(orbit, c);
            seen.(c) := true;
            if BottIsUpperTriangularAux() then
                AVLDelete( t, c );
            fi;
        fi;
    end;

    code := AVLIndexDelete( t, 1 );
    if code = fail then
        return fail;
    fi;
    orbit:= [ code ];

    seen := rec();
    seen.(code) := true;

    pos := 1;
    while pos<=Length(orbit) do

        BottDecodeMat( orbit[pos] );
        pos := pos+1;

        i := 0;
        while i<dim do
            i := i+1;
            j := i;
            while j<dim do
                j := j+1;
                BottSwapRowsAndCols(i,j);
                AddCode();
            od;
        od;

        i := 0;
        while i<dim do
            i := i+1;
            BottConditionalAddCol(i);
            AddCode();
        od;

        i := 0;
        while i<dim do
            i := i+1;
            j := i;
            while j < dim do
                j := j+1;
                if BottConditionalAddRow(i, j) then
                    AddCode();
                fi;
                if BottConditionalAddRow(j, i) then
                    AddCode();
                fi;
            od;
        od;
    od;

    return code;
end;

BottAVLNext := BottAVLNextPos;

BottAVLRepresentatives := function(t, dim)
    local r, l;
    l := [];
    r := BottAVLNext(t, dim);
    while r<>fail do
        Info( BottCalcInfo, 1, "[BOTT] Added representative with code: ", r );
        Info( BottCalcInfo, 2, "[BOTT] Left elements in the tree: ", Length(t));
        Add(l, r);
        r := BottAVLNext(t, dim);
    od;
    return l;
end;
