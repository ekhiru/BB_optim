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
                    INSTANCES="$INSTANCES LOP-synthetic_seed=${seed}_n=${n}_m=${m}_phi=${phi}"
                done
            done
        done
    done
    echo $INSTANCES
}

INSTANCES="\
pfsp/rec05.txt \
pfsp/rec13.txt \
pfsp/rec19.txt \
pfsp/rec31.txt \
"


###### For LOLIB instances
INSTANCES="$INSTANCES $(grep -v '#' lolib-instances.txt | tr '\n' ' ')"

###### Synthetic LOP instances
#INSTANCES=$(gen_lop_synthetic $INSTANCES)

OUTFILE="$BINDIR/results/m400-er0/tuning.csv"
echo "instance,seed,umm_m_ini,budgetMM,r_1,r_2,budget,eval_ranks,fitness" > $OUTFILE 

for instance in $INSTANCES; do
    RESULTS="$OUTDIR/tuning-er?/$instance"
    echo "Reading $RESULTS"
    cat $RESULTS/* >> $OUTFILE
done
echo "$OUTFILE created"
