from problem import Problem
import mallows_kendall as mk
import numpy as np
import re

def synthetic_LOP(n, m, phi):
  instance = np.zeros((n,n))
  s = np.asarray(mk.samplingMM(m, n, phi=phi, k=None))
  for i in range(n):
      for j in range(i+1,n):
          instance[i,j] = (s[:,i] < s[:,j]).sum() / m
          instance[j,i] = 1 - instance[i,j]
  return instance

# Linear Ordering Problem
class LOP(Problem):
  # Class attributes
  problem_name = "LOP"
  
  @classmethod
  def generate_synthetic(cls, n, m, phi, seed = None):
    if seed is None:
      seed = np.random.randint(1, 123456789, 1)[0]
    np.random.seed(seed)
    instance = synthetic_LOP(n, m, phi)
    best_sol = list(range(n))
    worst_sol = best_sol[::-1]
    return cls(n, instance, best_sol = best_sol, worst_sol = worst_sol,
               instance_name = f"LOP-synthetic,seed={seed},n={n},m={m},phi={phi}")
    
  @classmethod
  def read_instance(cls, filename):
    if "synthetic" in filename:
      seed, n, m, phi = re.search("seed=([0-9]+),n=([0-9]+),m=([0-9]+),phi=([^ ]+)", filename).group(1,2,3,4)
      return cls.generate_synthetic(int(n), int(m), float(phi), int(seed))
    else:
      # FIXME: How to read LOPLIB instances?
      pass
    
  # Methods
  def __init__(self, n, instance, best_sol = None, worst_sol = None, instance_name = None):
    self.n = n
    # FIXME: Can we find a better name than instance?
    self.instance = instance
    super().__init__(best_sol = best_sol, worst_sol = worst_sol, instance_name = instance_name)
        
  def fitness_nosave(self, x):
    # In case it is not numpy array.
    x = np.asarray(x, dtype=int)
    # Sum of the upper triangle. We have to minimize this.
    f = np.tril(self.instance[np.ix_(x, x)]).sum()
    return f

  
