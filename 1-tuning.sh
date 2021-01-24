#!/bin/bash
set -e
set -o pipefail

# Find our own location.
BINDIR=$(dirname "$(readlink -f "$(type -P $0 || echo $0)")")
OUTDIR="$HOME/scratch"

# This function launches one job $1 is the job name, the other arguments is the job to submit.
qsub_job() {
    PARALLEL_ENV=smp.pe
    # We would like to use $BASHPID here, but OS X version of bash does not
    # support it.
    JOBNAME=umm-$counter-r$run-$$
    OUTFILE="${RESULTS}/$JOBNAME"
    
    qsub -v PATH <<EOF
#!/bin/bash --login
#$ -N $JOBNAME
# -pe $PARALLEL_ENV $NB_PARALLEL_PROCESS 
#$ -l ivybridge
#$ -M manuel.lopez-ibanez@manchester.ac.uk
#$ -m ase
#      b     Mail is sent at the beginning of the job.
#      e     Mail is sent at the end of the job.
#      a     Mail is sent when the job is aborted or rescheduled.
#      s     Mail is sent when the job is suspended.
#
#$ -o $OUTDIR/${JOBNAME}.stdout
#$ -j y
#$ -cwd
module load apps/anaconda3
echo "running: $BINDIR/target-runner-umm.py umm $counter-r$run-$$ $run $instance --m_ini $umm_m_ini --budgetMM $budgetMM --rsl $r_1 --wml $r_2 --budget $budget --eval_ranks $eval_ranks"
echo -n "$instance,$run,$umm_m_ini,$budgetMM,$r_1,$r_2,$budget,$eval_ranks,$value" > $OUTFILE
$BINDIR/target-runner-umm.py umm $counter-r$run-$$ $run $instance --m_ini $umm_m_ini --budgetMM $budgetMM --rsl $r_1 --wml $r_2 --budget $budget --eval_ranks $eval_ranks >> $OUTFILE
EOF

}
launch_local() {
    JOBNAME=umm-$counter-r$run-$$
    OUTFILE=$RESULTS/$JOBNAME
    echo "running: $BINDIR/target-runner-umm.py umm $counter-r$run-$$ $run $instance --m_ini $umm_m_ini --budgetMM $budgetMM --rsl $r_1 --wml $r_2 --budget $budget --eval_ranks $eval_ranks"
    echo -n "$instance,$run,$umm_m_ini,$budgetMM,$r_1,$r_2,$budget,$eval_ranks,$value" > $OUTFILE
    $BINDIR/target-runner-umm.py umm $counter-r$run-$$ $run $instance --m_ini $umm_m_ini --budgetMM $budgetMM --rsl $r_1 --wml $r_2 --budget $budget --eval_ranks $eval_ranks >> $OUTFILE
}

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

nruns=10

LAUNCHER=qsub_job
#LAUNCHER=launch_local

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

budget=400

#eval_ranks=1
eval_ranks=0

r_1_values=$(seq 0.1 0.1 0.5)
r_2_values="$(seq 0.6 0.1 0.9) 0.99"
budgetMM=10
umm_m_ini=10

counter=0
for instance in $INSTANCES; do
    RESULTS="$OUTDIR/tuning-er${eval_ranks}/$instance"
    mkdir -p $RESULTS
    for r_1 in $r_1_values; do
	for r_2 in $r_2_values; do
	    counter=$((counter+1))
	    for run in $(seq 1 $nruns); do
		$LAUNCHER $counter $run $RESULTS $instance $umm_m_ini $budgetMM $r_1 $r_2 $budget $eval_ranks
	    done
	done
    done
done
