#!bin/bash

if [ -r umm.zip ]; then
    rm -f umm.zip
fi

cat << EOF | zip --exclude \*~ -r umm.zip -@
README.txt
0-setup.sh
1-tuning.sh
2-collect-tuning.sh
3-launch_exps.sh
4-analysis.ipynb
best_fitness_selected.csv
CEGO_2.4.0.tar.gz
cego.py
img
lolib-instances.txt
lop
lop.py
mallows_kendall.py
pfsp
pfsp.py
problem.py
qap.py
results
runner.py
target-runner-cego.py
target-runner-umm.py
umm.py
EOF


