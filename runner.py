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

def get_problem(instance_name):
    instance_name = instance_name.lower()
    if "qap" in instance_name:
        from qap import QAP
        return QAP
    elif "lop" in instance_name:
        from lop import LOP
        return LOP
    elif "pfsp" in instance_name:
        from pfsp import pfsp
        return PFSP
    raise ValueError("Unknown problem: " + instance_name)


def run_once(algo, instance, seed, **algo_params):
    timer = Timer()
    df = algo(instance, seed, **algo_params)
    if instance.best_fitness is not None and instance.worst_fitness is not None:
        df['Fitness'] = (df.Fitness - instance.best_fitness) / (instance.worst_fitness - instance.best_fitness)
    
    df['run_time'] = timer.elapsed()
    return df

# from cego import cego
# from lop import LOP
# from umm import uMM

# instance = LOP.read_instance("synthetic,seed=1,n=20,m=200,phi=0.5")

# #df = run_one(cego, instance, 1, budget = 10, m_ini = 5, budgetGA=100)
# df = run_one(uMM, instance, 1, budget = 100, m_ini = 5, budgetMM=10,
#              ratio_samples_learn = 0.22,
#              weight_mass_learn = .83)

