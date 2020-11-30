#!/bin/bash
set -e
set -o pipefail

# Find our own location.
BINDIR=$(dirname "$(readlink -f "$(type -P $0 || echo $0)")")

# This function launches one job $1 is the job name, the other arguments is the job to submit.
qsub_job() {
    PARALLEL_ENV=smp.pe
    # We would like to use $BASHPID here, but OS X version of bash does not
    # support it.
    JOBNAME=$1-$$
    shift 1
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
#$ -o ${JOBNAME}.stdout
#$ -j y
#$ -cwd
module load apps/anaconda3
echo "running: $@"
$@
EOF
}

launch_local() {
    echo $1
    shift
    echo "running: $@"
    $@
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
 qap/kra32.dat \
 qap/nug12.dat \
 qap/nug30.dat \
 qap/tho30.dat \
 pfsp/rec05.txt \
 pfsp/rec13.txt \
 pfsp/rec19.txt \
"
# pfsp/rec31.txt \

###### For LOLIB instances
#INSTANCES="$INSTANCES $(grep -v '#' lolib-instances.txt | tr '\n' ' ')"

###### Synthetic LOP instances
#INSTANCES=$(gen_lop_synthetic "")

budget=400

eval_ranks=0
#eval_ranks=1

cego_m_ini=10
budgetGA=3 # Actually, 10**budgetGA

r_1=0.1
r_2=0.9
budgetMM=10
umm_m_ini=10

counter=0
for instance in $INSTANCES; do
    counter=$((counter+1))
    RESULTS="results-er${eval_ranks}/$instance"
    mkdir -p $RESULTS
    for run in $(seq 1 $nruns); do
       	### Uncomment for running CEGO
	$LAUNCHER cego-$counter-r$run ./target-runner-cego.py cego $counter-r$run-$$ $run $instance --m_ini $cego_m_ini --budgetGA $budgetGA --budget $budget --eval_ranks $eval_ranks --output $RESULTS/cego-r$run

	### Uncomment for running UMM
	$LAUNCHER umm-$counter-r$run ./target-runner-umm.py umm $counter-r$run-$$ $run $instance --m_ini $umm_m_ini --budgetMM $budgetMM --rsl $r_1 --wml $r_2 --budget $budget  --eval_ranks $eval_ranks --output $RESULTS/umm-r$run
        
    done
done

RERUN_params=0
if [ $RERUN_params -eq 0 ]; then
    exit 0
fi

counter=0
for instance in $INSTANCES; do
    for r_1 in $(seq 0.1 0.1 0.5); do
	for r_2 in "$(seq 0.6 0.1 0.9) 0.99"; do
	    counter=$((counter+1))
	    RESULTS="results-r1r2/$instance"
	    mkdir -p $RESULTS
	    for run in $(seq 1 $nruns); do
		### Uncomment for running UMM
		$LAUNCHER umm-$counter-r$run ./target-runner-umm.py umm $counter-r$run-$$ $run $instance --m_ini $umm_m_ini --budgetMM $budgetMM --rsl $r_1 --wml $r_2 --budget $budget --eval_ranks $eval_ranks --output $RESULTS/umm-r$run
	    done
	done
    done
done
