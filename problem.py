import numpy as np
from rpy2.robjects import FloatVector
import rpy2.rinterface as ri

from mallows_kendall import kendallTau

class Problem:
    def __init__(self, best_sol = None, worst_sol = None, instance_name = None, best_fitness=None):
        self.best_sol = best_sol
        self.worst_sol = worst_sol
        self.instance_name = instance_name
        if best_sol is None:
            self.best_fitness = best_fitness
        else:
            self.best_fitness = self.fitness_nosave(best_sol)
        if worst_sol is None:
            self.worst_fitness = None
        else:
            self.worst_fitness = self.fitness_nosave(worst_sol)
        self.reset()

    def reset(self):
        self.evaluations = []
        self.solutions = []

    def fitness_nosave(self, x):
        raise NotImplementedError("virtual method")

    def fitness(self, x):
      f = self.fitness_nosave(x)
      self.solutions.append(x)
      self.evaluations.append(f)
      return f

    # FIXME: Is this distance correct for all problems?
    def distance_to_best(self, perm):
        if self.best_sol is None:
            return np.nan
        return kendallTau(perm, self.best_sol) / (self.n * (self.n - 1) * 0.5)

    # Returns a closure function that can be called from R.
    # WARNING: this function minimizes for CEGO
    def make_r_fitness(self):
        @ri.rternalize
        def r_fitness(x):
            xpy = np.asarray(x) - 1 # R vectors are 1-indexed
            y = self.fitness(xpy)
            return FloatVector(np.asarray(y))
        return r_fitness
