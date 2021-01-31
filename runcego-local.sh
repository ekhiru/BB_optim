#!/bin/bash
parallel --tag R --slave --no-save --no-restore -f ./runCEGO.R --args {2} {1} ::: 0 1 ::: $(seq 1 10)
