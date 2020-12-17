#!/bin/bash --login
#$ -cwd
#$ -t 1-60
#$ -l ivybridge
#$ -M manuel.lopez-ibanez@manchester.ac.uk
#$ -m ase
#      b     Mail is sent at the beginning of the job.
#      e     Mail is sent at the end of the job.
#      a     Mail is sent when the job is aborted or rescheduled.
#      s     Mail is sent when the job is suspended.
#
#$ -j y
# These should be used with caution (see below for explanation)
if [[ $SGE_TASK_ID -gt 30 ]]; then
    SEED=$((SGE_TASK_ID-30))
    EVAL_RANKS=1
else
    SEED=$SGE_TASK_ID
    EVAL_RANKS=0
fi
R --slave --no-save --no-restore -f ./runCEGO.R --args $SEED $EVAL_RANKS

