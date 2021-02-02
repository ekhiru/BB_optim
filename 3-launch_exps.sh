#!/bin/bash
set -e
set -o pipefail

# Find our own location.
BINDIR=$(dirname "$(readlink -f "$(type -P $0 || echo $0)")")
OUTDIR="$HOME/scratch-new"

# This function launches one job $1 is the job name, the other arguments is the job to submit.
qsub_job() {
    PARALLEL_ENV=smp.pe
    # We would like to use $BASHPID here, but OS X version of bash does not
    # support it.
    ALGO=$1
    OUTPUT=$2
    shift 2
    JOBNAME=${ALGO}-$counter-$$
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
echo "running: ./target-runner-${ALGO}.py $ALGO $counter-$$-r\$run \$run $@ --output ${OUTPUT}-r\$run"
./target-runner-${ALGO}.py $ALGO $counter-$$-r\$run \$run $@ --output "${OUTPUT}-r\$run"
EOF
}

# FIXME: Not working right now
launch_local() {
    ALGO=$1
    OUTPUT=$2
    shift 2
    for run in $(seq 1 $nruns); do
	echo "running: ./target-runner-${ALGO}.py $ALGO $counter-$$-r$run $run $@ --output ${OUTPUT}-r$run"
	./target-runner-${ALGO}.py $ALGO $counter-$$-r$run $run $@ --output "${OUTPUT}-r$run"
    done
}

# Generate LOP synthetic
gen_lop_synthetic() {
    instances=""
    LOP_n=20
    LOP_m=200
    LOP_seed=123456
    LOP_phi="0.5 0.7 0.9"

    for n in $LOP_n; do
        for m in $LOP_m; do
            for seed in $LOP_seed; do
                for phi in $LOP_phi; do
                    instances="$instances LOP-synthetic,seed=${seed},n=${n},m=${m},phi=${phi}"
                done
            done
        done
    done
    echo $instances
}

nruns=10

LAUNCHER=qsub_job
#LAUNCHER=launch_local

## For QAP, PFSP instances
INSTANCES="\
  qap/kra30a.dat \
  qap/kra30b.dat \
 qap/kra32.dat \
 qap/nug12.dat \
 qap/nug30.dat \
 qap/tho30.dat \
 pfsp/rec05.txt \
 pfsp/rec13.txt \
pfsp/rec19.txt \
pfsp/rec31.txt \
"

INSTANCES="\
pfsp/rec05.txt \
pfsp/rec13.txt \
pfsp/rec19.txt \
pfsp/rec31.txt \
"


###### For LOLIB instances
INSTANCES="$INSTANCES $(grep -v '#' lolib-instances.txt | tr '\n' ' ')"
#INSTANCES="$(grep -v '#' lolib-instances.txt | tr '\n' ' ')"

###### Synthetic LOP instances
#INSTANCES="$INSTANCES $(gen_lop_synthetic)"

budget=400

#eval_ranks=1
eval_ranks="0 1"
eval_ranks=0

cego_m_ini=10
# Actually, 10**budgetGA
budgetGA=4

r_1=0.1
r_2=0.9
#budgetMM=10
budgetMM=1
umm_m_ini=10

counter=0
for er in $eval_ranks; do
    for instance in $INSTANCES; do
	counter=$((counter+1))
	RESULTS="$OUTDIR/results/m${budget}-er${eval_ranks}/$instance"
	mkdir -p "$RESULTS"
	$LAUNCHER umm "${RESULTS}/umm-b${budgetMM}" $instance --m_ini $umm_m_ini --budgetMM $budgetMM --rsl $r_1 --wml $r_2 --budget $budget --eval_ranks $er
	#$LAUNCHER cego "${RESULTS}/cego" $instance --m_ini $cego_m_ini --budgetGA $budgetGA --budget $budget --eval_ranks $er
    done
done


