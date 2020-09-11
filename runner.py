import time
class Timer:
    def __init__(self):
        self.now = time.perf_counter()

    def elapsed(self):
        old_now = self.now
        self.now = time.perf_counter()
        return self.now - old_now

    def reset(self):
        self.elapsed()

import pandas as pd
import numpy as np

def get_instance(instance_name):
    instance_name_lower = instance_name.lower()
    if "qap" in instance_name_lower:
        from qap import QAP
        return QAP.read_instance(instance_name, opt_filename = "qap/optimal.txt")
    elif "lop" in instance_name_lower:
        from lop import LOP
        return LOP.read_instance(instance_name)
    elif "pfsp" in instance_name_lower:
        from pfsp import PFSP
        return PFSP.read_instance(instance_name)
    
    raise ValueError("Unknown problem: " + instance_name)


def run_once(algo_name, instance_name, seed, out_filename = None,
             **algo_params):
    if algo_name == "uMM":
        from umm import uMM
        algo = uMM
    elif algo_name == "CEGO":
        from cego import cego
        algo = cego
    else:
        raise ValueError("Unknown algo: " + algo_name)
    
    instance = get_instance(instance_name)
    timer = Timer()
    df = algo(instance, seed, **algo_params)
    if instance.best_fitness is not None and instance.worst_fitness is not None:
        df['Fitness'] = (df.Fitness - instance.best_fitness) / (instance.worst_fitness - instance.best_fitness)
    df['Function evaluations'] = np.arange(1, len(df['Fitness'])+1)
    df['run_time'] = timer.elapsed()
    df['Problem'] = instance.problem_name
    df['instance'] = instance.instance_name
    df['Solver'] = algo_name
    if out_filename is not None:
        df.to_csv(out_filename + '.csv.gz', index=False)
        df.to_pickle(out_filename + '.pkl.gz')
    return df

