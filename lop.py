import numpy as np
import mallows_kendall as mk

from rpy2.robjects import FloatVector
import rpy2.rinterface as ri

def synthetic_LOP(n, m, phi):
  instance = np.zeros((n,n))
  s = np.asarray(mk.samplingMM(m,n, phi=phi, k=None))
  for i in range(n):
      for j in range(i+1,n):
          instance[i,j] = (s[:,i]< s[:,j]).sum()
          instance[j,i] = m - instance[i,j]
  return instance

def get_fitness(perm, instance):
  # In case it is not numpy array
  perm = np.asarray(perm)
  n = len(perm) #sum in the upper triangle. we have to maximize this
  sol = 0
  for i in range(n):
    sol += instance[perm[i], perm[i:n]].sum()
  return sol

# Linear Ordering Problem
class LOP:
  # Class attributes
  problem_name = "LOP"

  @classmethod
  def generate_synthetic(cls, n, m, phi):
    return cls(n, synthetic_LOP(n, m, phi))
    
  # Methods
  def __init__(self, n, instance):
    self.n = n
    self.instance = instance
    self.evaluations = []
    self.solutions = []

  def fitness_nosave(self, perm):
      return get_fitness(perm, self.instance)

  # an alias of the above for compatibility
  def get_fitness(self, perm):
    return self.fitness_nosave(perm)
      
  # Minimized
  def fitness(self, perm):
      x = get_fitness(perm, self.instance)
      self.evaluations.append(x)
      self.solutions.append(perm)
      return x

  # Returns a closure function that can be called from R.
  # WARNING: this function minimizes for CEGO
  def make_r_fitness(self):
      @ri.rternalize
      def r_fitness(x):
          y = self.fitness(x)
          return FloatVector(-np.asarray(y))
      return r_fitness
