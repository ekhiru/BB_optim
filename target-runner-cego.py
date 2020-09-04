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
import pandas as pd
from cego import cego
import runner

from argparse import ArgumentParser,RawDescriptionHelpFormatter,_StoreTrueAction,ArgumentDefaultsHelpFormatter,Action
parser = ArgumentParser(description = "CEGO")
parser.add_argument('configuration_id', type=str, help='configuration_id')
parser.add_argument('instance_id', type=str, help='instance_id')
parser.add_argument('seed', type=int, help='random seed')
parser.add_argument('instance_name', type=str, help='instance name')

# Parameters for the target algorithm
parser.add_argument('--m_ini', type=int, default=0, help='m_ini')
parser.add_argument('--budget', type=int, default=400, help='budget')
parser.add_argument('--budgetGA', type=int, default=0, help='budgetGA')
parser.add_argument("--output", type=str, default=None, help="output file")

args = parser.parse_args()

problem = runner.get_problem(args.instance_name)
instance = problem.read_instance(args.instance_name)

budget = args.budget
assert budget > 2 * args.m_ini

budgetGA = 10**args.budgetGA

stdout = sys.stdout
outfilename = f'c{args.configuration_id}-{args.instance_id}-{args.seed}.stdout' 
with open(outfilename, 'w') as sys.stdout:
    df = runner.run_once(cego, instance, args.seed, budget = budget, m_ini = args.m_ini, budgetGA = budgetGA)
    if args.output is not None:
        df['Problem'] = instance.problem_name
        df['instance'] = args.instance_name
        df['Solver'] = "CEGO"
        df.to_csv(args.output + '.csv', index=False)
        df.to_pickle(args.output + '.pkl')
        
sys.stdout = stdout
print(df["Fitness"].min())
# remove tmp file.
os.remove(outfilename)
