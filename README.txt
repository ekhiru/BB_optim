<!-- Best viewed as markdown                        -*- mode: markdown -*- -->

Reproducibility steps for the paper
===================================

Setup
-----

Read first, then run the `0-setup.sh` script.


Re-run parameter analysis
-------------------------

The script `1-tuning.sh` creates the data for the parameter analysis. You may
need to adapt the script for your computing system.

The data is collected into one CSV file by the script `2-collect-tuning.sh`.


Re-run experiments
------------------

The script `3-launch-exps.sh` launches the main experiments. 

You have to select the parameter settings, instances and whether to run CEGO
and/or UMM. You may need to adapt the script for your computing system.


How to reproduce the plots and tables without re-running
---------------------------------------------------------

In notebook `4-analysis.ipynb`, evaluate the first cell (imports),
then evaluate the cells in the sections that interest you.

