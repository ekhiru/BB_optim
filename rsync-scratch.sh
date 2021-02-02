#!/bin/bash
rsync -av $HOME/umm $HOME/scratch/ --exclude-from=- <<"EOF"
.git
results*
old
pickles
paper
*~
img/
irace/
__pycache__
rsync-scratch.sh
EOF




