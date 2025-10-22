#!/bin/bash

ok=0
fail=0

start=$(date +%s%N)   # nanoseconds
for exe in $(find ./ -maxdepth 1 -type f -perm /111 -regex './test-[^\.]*'); do
    ./$exe
    status=$?
    if [ $status -eq 0 ]; then
        ((ok++))
    else
        ((fail++))
    fi
    echo ""
done
end=$(date +%s%N)     # nanoseconds
elapsed=$((end - start))
elapsed_s=$(awk "BEGIN {printf \"%.3f\", $elapsed / 1000000000}")

echo   "=== Summary: ==="
echo   "    passed: $ok test(s)"
echo   "    failed: $fail test(s)"
echo   "    time:   $elapsed_s s"