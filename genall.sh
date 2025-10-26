#!/bin/bash

NJOBS=$(nproc --all --ignore=2)

OUTDIR=.
if test -d $1; then
    OUTDIR=$1
fi;

echo "Setting output directory to '$OUTDIR'"

for ((d=1; d<=8; d++)); do
    # echo "Starting dimension $d."
    dag=$OUTDIR/0$d-dag-all.d6
    geng -q $d | directg -qa > $dag
    bott=$OUTDIR/0$d-rbm-all.d6
    ./minimalf -i $dag -o $bott
    orientable=$OUTDIR/0$d-rbm-orientable.d6
    cat $bott | ./orientedf > $orientable
    spinc=$OUTDIR/0$d-rbm-spinc.d6
    cat $orientable | ./upperg | ./spincf > $spinc
    spin=$OUTDIR/0$d-rbm-spin.d6
    cat $spinc | ./spinf > $spin
    echo "Summary of dimension $d:"
    echo "  RBM in total:   $(cat $bott | wc -l)"
    echo "  Orientable RBM: $(cat $orientable | wc -l)"
    echo "  Spinc RBM:      $(cat $spinc | wc -l)"
    echo "  Spin RBM:       $(cat $spin | wc -l)"
done

graphs=$OUTDIR/09.g6
geng -q 9 > $graphs
echo "[INFO] Generating DAGs with 9 vertices ..." >&2
# dag=$OUTDIR/09-dag-all.xz
# pv $graphs | directg -qa | xz -T $NJOBS -c > $dag

# dag=$OUTDIR/09-dag-all.d6
# pv $graphs | directg -qa > $dag

split -l 2500 $graphs $OUTDIR/9_
files=$(find $OUTDIR -regex '.*/9_[a-z]+')
parallel -j $NJOBS --lb 'directg -a {} {}.out && rm {}' ::: $files
dag=$OUTDIR/09-dag-all.d6
cat $OUTDIR/9_*.out > $dag
rm -f $OUTDIR/9_*.out
bott=$OUTDIR/09-rbm-all.d6
echo "[INFO] Generating RBMs of dimension 9 ..." >&2
# xzcat $dag| ./minimalf -o $bott -l 100M -n 1023 -vv -j $NJOBS
./minimalf -i $dag -o $bott -l 100M -n 1023 -vv -j $NJOBS
orientable=$OUTDIR/09-rbm-orientable.d6
cat $bott | ./orientedf > $orientable
spinc=$OUTDIR/09-rbm-spinc.d6
cat $orientable | ./upperg | ./spincf > $spinc
spin=$OUTDIR/09-rbm-spin.d6
cat $spinc | ./spinf > $spin
printf "\
Summary of dimension 9:\n\
  RBM in total:   $(cat $bott | wc -l)\n\
  Orientable RBM: $(cat $orientable | wc -l)\n\
  Spinc RBM:      $(cat $spinc | wc -l)\n\
  Spin RBM:       $(cat $spin | wc -l)\n"

### This approach is possible, however may be done faster
# graphs=$OUTDIR/10.g6
# geng -q 10 > $graphs
# split -l 25000 $graphs $OUTDIR/10_
# files=$(find $OUTDIR -regex '.*/10_[a-z]+')
### the parallel job below took about 10 hours and 50 minutes (i7-10700)
# parallel -j $((NJOBS * 2 / 3)) --lb 'directg -a {} | ./orientedf > {}.out && rm {}' ::: $files
dag=$OUTDIR/10-dag-orientable.d6
### the orientedg, which did the same as parallel above, took 1 hour and 18 minutes
echo "[INFO] Generating DAGs with 10 vertices all of even degrees ..." >&2
./orientedg -d 10 -o $dag -p
orientable=$OUTDIR/10-rbm-orientable.d6
echo "[INFO] Generating orientable RBMs of dimension 10 ..." >&2
./minimalf -i $dag -o $orientable -l 100M -n 1023 -vv -j $NJOBS
spinc=$OUTDIR/10-rbm-spinc.d6
cat $orientable | ./upperg | ./spincf > $spinc
spin=$OUTDIR/09-rbm-spin.d6
cat $spinc | ./spinf > $spin
printf "Summary of dimension 10:\n\
  Orientable RBM: $(cat $orientable | wc -l)\n\
  Spinc RBM:      $(cat $spinc | wc -l)\n\
  Spin RBM:       $(cat $spin | wc -l)\n"

dag=$OUTDIR/11-dag-spinc.d6
echo "[INFO] Generating \"spinc\" DAGs with 11 vertices ..." >&2
./backtrack -d 11 -s 9 -j $NJOBS -o $dag -p -n -v
spinc=$OUTDIR/11-rbm-spinc.d6
echo "[INFO] Generating spinc RBMs of dimension 11 ..." >&2
./minimalf -i $dag -o $spinc -l 100M -n 1023 -vv -j 22
spin=$OUTDIR/11-rbm-spin.d6
cat $spinc | ./spinf > $spin
printf "Summary of dimension 11:\n\
  Spinc RBM:      $(cat $spinc | wc -l)\n\
  Spin RBM:       $(cat $spin | wc -l)\n"