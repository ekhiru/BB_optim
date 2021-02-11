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
    JOBNAME=umm-$counter-$$
    qsub -v PATH <<EOF
#!/bin/bash --login
#$ -t 1-$nruns
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
run=\$SGE_TASK_ID
OUTFILE="${RESULTS}/$JOBNAME-r\$run"
echo "running: $BINDIR/target-runner-umm.py umm $counter-$$-r\$run \$run $@ > \$OUTFILE"
echo -n "\"$instance\",\$run,$umm_m_ini,$beta,$r_1,$r_2,$budget,$eval_ranks," > \$OUTFILE
$BINDIR/target-runner-umm.py umm $counter-$$-r\$run \$run $@ >> \$OUTFILE
EOF

}
launch_local() {
    for run in $(seq 1 $nruns); do
	JOBNAME=umm-$counter-$$
	OUTFILE=$RESULTS/${JOBNAME}-r$run
	echo "running: $BINDIR/target-runner-umm.py umm $counter-$$-r$run $run $@ > $OUTFILE"
	echo -n "\"$instance\",$run,$umm_m_ini,$beta,$r_1,$r_2,$budget,$eval_ranks," > $OUTFILE
	$BINDIR/target-runner-umm.py umm $counter-$$-r$run $run $@ >> $OUTFILE
    done
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
                    INSTANCES="$INSTANCES LOP-synthetic_seed=${seed}_n=${n}_m=${m}_phi=${phi}"
                done
            done
        done
    done
    echo $INSTANCES
}

nruns=10

LAUNCHER=qsub_job
#LAUNCHER=launch_local

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

budget=400

#eval_ranks=1
eval_ranks=0

r_1_values=$(seq 0.1 0.1 0.5)
r_2_values="$(seq 0.6 0.1 0.9) 0.99"
beta=1
umm_m_ini=10

counter=0
for instance in $INSTANCES; do
    RESULTS="$OUTDIR/tuning-er${eval_ranks}/$instance"
    mkdir -p $RESULTS
    for r_1 in $r_1_values; do
	for r_2 in $r_2_values; do
	    counter=$((counter+1))
	    $LAUNCHER $instance --m_ini $umm_m_ini --budgetMM $beta --rsl $r_1 --wml $r_2 --budget $budget --eval_ranks $eval_ranks
	done
    done
done
