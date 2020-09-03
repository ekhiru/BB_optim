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
    JOBNAME=$$-$1
    shift 1
    qsub -v PATH <<EOF
#!/bin/bash --login
#$ -N $JOBNAME
#$ -pe $PARALLEL_ENV $NB_PARALLEL_PROCESS 
#$ -M manuel.lopez-ibanez@manchester.ac.uk
#$ -m ase
#      b     Mail is sent at the beginning of the job.
#      e     Mail is sent at the end of the job.
#      a     Mail is sent when the job is aborted or rescheduled.
#      s     Mail is sent when the job is suspended.
#
#$ -o ${JOBNAME}.stdout
#$ -j
#$ -cwd
module load apps/R/3.5.2
module load apps/anaconda3
$@
EOF
}

n=20
budgetGA=10
budgetMM=10
m_max=400
phi_instances="0.5 0.7 0.9"
nruns = 9 # we run +1

LAUNCHER=qsub_job

for repe in $(seq 0 $nruns); do
  for phi_instance in $phi_instances; do
      ratio_samples_learn=0.15
      weight_mass_learn=0.84
      budgetMM=15
      m_ini=23
      $LAUNCHER python3 cego_lop.py $n $repe $phi_instance $budgetGA $budgetMM $ratio_samples_learn $weight_mass_learn $SLURM_JOB_ID $m_max $m_ini

      ratio_samples_learn=0.23
      weight_mass_learn=0.88
      budgetMM=16
      m_ini=38
      #sbatch --export=n=$n,repe=$repe,phi_instance=$phi_instance,budgetGA=$budgetGA,m_max=$m_max,budgetMM=$budgetMM,ratio_samples_learn=$ratio_samples_learn,weight_mass_learn=$weight_mass_learn,m_ini=$m_ini run_one.sl
      $LAUNCHER python3 cego_lop.py $n $repe $phi_instance $budgetGA $budgetMM $ratio_samples_learn $weight_mass_learn $SLURM_JOB_ID $m_max $m_ini

      ratio_samples_learn=0.22
      weight_mass_learn=0.83
      budgetMM=14
      m_ini=32
      #sbatch --export=n=$n,repe=$repe,phi_instance=$phi_instance,budgetGA=$budgetGA,m_max=$m_max,budgetMM=$budgetMM,ratio_samples_learn=$ratio_samples_learn,weight_mass_learn=$weight_mass_learn,m_ini=$m_ini run_one.sl
      $LAUNCHER python3 cego_lop.py $n $repe $phi_instance $budgetGA $budgetMM $ratio_samples_learn $weight_mass_learn $SLURM_JOB_ID $m_max $m_ini

      ratio_samples_learn=0.17
      weight_mass_learn=0.83
      budgetMM=17
      m_ini=32
      #sbatch --export=n=$n,repe=$repe,phi_instance=$phi_instance,budgetGA=$budgetGA,m_max=$m_max,budgetMM=$budgetMM,ratio_samples_learn=$ratio_samples_learn,weight_mass_learn=$weight_mass_learn,m_ini=$m_ini run_one.sl
      $LAUNCHER python3 cego_lop.py $n $repe $phi_instance $budgetGA $budgetMM $ratio_samples_learn $weight_mass_learn $SLURM_JOB_ID $m_max $m_ini
    done
done

