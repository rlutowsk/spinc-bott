AddCode := function(seen, rval, orbit)
    local code;
        code := BottEncodeAux();
        if not IsBound(seen.(code)) then
            Add(orbit, code);
            seen.(code) := BottIsUpperTriangularAux();
            if seen.(code) then
                rval.(code) := true;
            fi;
        fi;
    end;


BottOrbit := function(dim)
    local orbit, x, i, j, code, seen, rval, first;

    AddCode := function()
        code := BottEncodeAux();
        if not IsBound(seen.(code)) then
            Add(orbit, code);
            seen.(code) := BottIsUpperTriangularAux();
            if seen.(code) then
                rval.(code) := true;
            fi;
        fi;
    end;

    code := BottEncodeMat();
    first:= code;
    orbit := [ code ];
    seen := rec();
    seen.(code) := BottIsUpperTriangularMat();
    rval := rec();
    if seen.(code) then
        rval.(code) := true;
    fi;

    while orbit<>[] do
        BottDecodeMat(Remove(orbit));
        for i in [1..dim] do
            for j in [i+1..dim] do
                BottSwapRowsAndCols(i, j);
                AddCode();
            od;
        od;
        for i in [1..dim] do
            BottConditionalAddCol(i);
            AddCode();
            for j in [i+1..dim] do
                BottConditionalAddRow(i, j);
                AddCode();
                BottConditionalAddRow(j, i);
                AddCode();
            od;
        od;
    od;
    BottDecodeMat( first );
    return rval;
end;
