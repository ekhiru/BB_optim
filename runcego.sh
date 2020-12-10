#!/bin/bash --login
#$ -cwd
#$ -t 1-30
#$ -l ivybridge
#$ -M manuel.lopez-ibanez@manchester.ac.uk
#$ -m ase
#      b     Mail is sent at the beginning of the job.
#      e     Mail is sent at the end of the job.
#      a     Mail is sent when the job is aborted or rescheduled.
#      s     Mail is sent when the job is suspended.
#
#$ -o cegorun.stdout
#$ -j y
R --slave --no-save --no-restore -f ./runCEGO.R --args $SGE_TASK_ID
