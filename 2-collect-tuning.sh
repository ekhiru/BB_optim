#!/bin/bash
set -e
set -o pipefail

# Find our own location.
BINDIR=$(dirname "$(readlink -f "$(type -P $0 || echo $0)")")
OUTDIR="$HOME/scratch"

# Generate LOP synthetic
gen_lop_synthetic() {
    INSTANCES=$@
    LOP_n=20
    LOP_m=200
    LOP_seed=123456
    LOP_phi="0.5 0.7 0.9"

    for n in $LOP_n; do
        for m in $LOP_m; do
            for seed in $LOP_seed; do
                for phi in $LOP_phi; do
                    INSTANCES="$INSTANCES LOP-synthetic,seed=${seed},n=${n},m=${m},phi=${phi}"
                done
            done
        done
    done
    echo $INSTANCES
}

## For QAP, PFSP instances
INSTANCES="\
qap/kra30a.dat \
qap/kra30b.dat \
qap/nug30.dat \
qap/tho30.dat \
pfsp/rec13.txt \
pfsp/rec19.txt \
lop/RandB/N-p40-02 \
lop/IO/N-t59d11xx \
lop/SGB/N-sgb75.02 \
lop/xLOLIB/N-be75eec_150 \
"

###### Synthetic LOP instances
INSTANCES=$(gen_lop_synthetic $INSTANCES)

OUTFILE="$BINDIR/results-er0/tuning.csv"
echo "instance,seed,umm_m_ini,budgetMM,r_1,r_2,budget,eval_ranks,fitness" > $OUTFILE 

for instance in $INSTANCES; do
    RESULTS="$OUTDIR/tuning-er?/$instance"
    echo "Reading $RESULTS"
    cat $RESULTS/* >> $OUTFILE
done
echo "$OUTFILE created"
