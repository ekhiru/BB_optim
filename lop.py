import numpy as np
import itertools as it
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

def u_phi(sample, s0, ws):
    m , n = np.asarray(sample).shape
    #if s0 is None: s0 = np.argsort(np.argsort(rankings.sum(axis=0))) #borda
    dist_avg = np.asarray([mk.kendallTau(perm, s0) for perm in sample]*ws).sum()/ws.sum() #np.mean(np.array([kendallTau(s0, perm) for perm in rankings]))
    try:
        theta = optimize.newton(mk.mle_theta_mm_f, 0.01, fprime=mk.mle_theta_mm_fdev, args=(n, dist_avg), tol=1.48e-08, maxiter=500, fprime2=None)
    except:
        #if dist_avg == 0.0: return s0, np.exp(-5)#=phi
        print("error. fit_mm. dist_avg=",dist_avg, dist_avg == 0.0)
        print(s0)
        raise
    if theta < 0:
        theta = 0.001
    return np.exp(-theta) # theta = - np.log(phi)

def get_fitness(perm, instance, problem):
    sol = 0
    n = len(perm)#sum inthe upper triangle. we have to maximize this
    inverse = np.argsort(perm)
    for i in range(n):
        for j in range(i,n):
            sol += instance[int(inverse[i]),int(inverse[j])]
    return sol

def uborda(sample, ws):
    mul = (sample * ws[:, None]).sum(axis=0)
    return np.argsort(np.argsort(mul))

def paint(df):
    for problem in df.problem.drop_duplicates():
      sns.set_style("whitegrid")
      color_variable = 'rho'
      y_variables = ['fitnesses','phi_estim','phi_sample','dist']
      for y_variable in y_variables:
        plt.figure(figsize=(15,5))
        palette = sns.color_palette("husl", len(df[color_variable].drop_duplicates()))
        sns.lineplot(x='m',hue=color_variable, y=y_variable, data=df[df.problem==problem], palette=palette)
        plt.show()

def get_expected_distance(iter_ratio, ini_dist):
    return ini_dist - ini_dist * iter_ratio
#get_expected_distance(0, 4)


class LOP:
  def __init__(self, n, instance):
    self.n = n
    self.instance = instance
    self.evaluations = []
    self.solutions = []

  def fitness_nosave(self, perm):
      return get_fitness(perm, self.instance, "LOP")
      
  # Minimized
  def fitness(self, perm):
      x = get_fitness(perm, self.instance, "LOP")
      self.evaluations.append(x)
      self.solutions.append(perm)
      return x

  # Returns a closure function that can be called from R.
  def make_r_fitness(self):
      @ri.rternalize
      def r_fitness(x):
          y = self.fitness(x)
          return FloatVector(np.asarray(y))
      return r_fitness
