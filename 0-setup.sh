#!/bin/sh

# Install the R package `CEGO` v2.4.0. Original is available at:
#   https://CRAN.R-project.org/package=CEGO
#
# We used R 3.4.4 for our experiments but later versions should also work.
R CMD INSTALL CEGO_2.4.0.tar.gz

# Install rpy2, we used rpy2-3.3.5
pip install --force-reinstall --update --user rpy2

