#!/usr/bin/env python3
###############################################################################
# This script is the command that is executed every run.
# Check the examples in examples/
#
# This script is run in the execution directory (execDir, --exec-dir).
#
# PARAMETERS:
# argv[1] is the candidate configuration number
# argv[2] is the instance ID
# argv[3] is the seed
# argv[4] is the instance name
# The rest (argv[5:]) are parameters to the run
#
# RETURN VALUE:
# This script should print one numerical value: the cost that must be minimized.
# Exit with 0 if no error, with 1 in case of error
###############################################################################
import sys
import os
import numpy as np
import pandas as pd
import mallows_kendall as mk
from lop import LOP
from cego_lop import solve_one_umm
import re


from argparse import ArgumentParser,RawDescriptionHelpFormatter,_StoreTrueAction,ArgumentDefaultsHelpFormatter,Action
parser = ArgumentParser(description = "CEGO")
parser.add_argument('configuration_id', type=int, help='configuration_id')
parser.add_argument('instance_id', type=int, help='instance_id')
parser.add_argument('algo_seed', type=int, help='random seed')
# These are part of the instance definition
parser.add_argument('-seed', type=int, help='inst_seed')
parser.add_argument('-m', type=int, help='inst_m')
parser.add_argument('-n', type=int, help='inst_n')
parser.add_argument('-phi', type=float, help='inst_n')
# Parameters for the target algorithm
parser.add_argument('--rsl', type=float, default=0, help='rsl')
parser.add_argument('--wml', type=float, default=0, help='wml')
parser.add_argument('--m_ini', type=int, default=0, help='m_ini')
parser.add_argument('--budgetMM', type=int, default=0, help='budgetMM')

args = parser.parse_args()
# inst_seed, inst_m, inst_n, phi = re.search("seed=([0-9]+)\s+m=([0-9]+)\s+n=([0-9]+)\s+phi=([^ ]+)", args.instance).group(1,2,3,4)
# inst_seed = int(inst_seed)
# inst_m = int(inst_m)
# inst_n = int(inst_n)
# phi = float(phi)

np.random.seed(args.seed)
instance = LOP.generate_synthetic(args.n, args.m, args.phi)

budget = 400

stdout = sys.stdout
outfilename = f'c{args.configuration_id}-{args.instance_id}-{args.algo_seed}.stdout' 
with open(outfilename, 'w') as sys.stdout:
    out = solve_one_umm(instance, budget, args.algo_seed, args.m_ini, instance.best_sol, instance.worst_sol, args.budgetMM, args.rsl, args.wml)
    
sys.stdout = stdout
print(np.min(out["Fitness"]))
# remove tmp file.
os.remove(outfilename)
