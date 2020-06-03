from rpy2.robjects.packages import STAP
import numpy as np
import itertools as it
import mallows_kendall as mk
import os
os.environ['RPY2_CFFI_MODE'] = "API" # bug in cffi 1.13.0 https://bitbucket.org/rpy2/rpy2/issues/591/runtimeerror-found-a-situation-in-which-we

from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
from rpy2.robjects import numpy2ri
from rpy2.robjects import FloatVector
numpy2ri.activate()
import rpy2.rinterface as ri


# funcion de distancia entre permutaciones
# MANUEL: There is a tentative implementation for scipy here: https://github.com/scipy/scipy/pull/7205/files
# but it looks quite different.
def kendallTau(A, B):
    n = len(A)
    pairs = it.combinations(range(n), 2)
    distance = 0
    for x, y in pairs:
        a = A[x] - A[y]
        try:
            b = B[x] - B[y]# if discordant (different signs)
        # MANUEL: when can such error happen?
        except:
            print("ERROR kendallTau, check b",A, B, x, y)
        if (a * b < 0):
            distance += 1
    return distance


@ri.rternalize
def r_kendallTau(a,b): return kendallTau(a,b)

def synthetic_LOP(n, m, phi):
  instance = np.zeros((n,n))
  s = np.asarray(mk.samplingMM(m,n, phi=phi, k=None))
  for i in range(n):
      for j in range(i+1,n):
          instance[i,j] = (s[:,i]< s[:,j]).sum()
          instance[j,i] = m - instance[i,j]
  return instance

def u_phi(sample,s0, ws):
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

def uborda(sample,ws):
    mul = (sample*ws[:, None]).sum(axis=0)
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
  def __init__(self, n, m, phi):
    self.n = n
    self.instance = synthetic_LOP(n, m, phi)

  # Minimized
  def fitness(self, perm):
      return get_fitness(perm, self.instance, "LOP")

  # Returns a closure function that can be called from R.
  def make_r_fitness(self):
      @ri.rternalize
      def r_fitness(x):
          y = self.fitness(x)
          return FloatVector(np.asarray(y))
      return r_fitness


def runCEGO(lop, mi, budget):
    rstring = """
    library(CEGO)
    my_cego <- function(fun, dist, n, mi = 5, budget = 15, seed = 0)
    {
    set.seed(seed)
    # mutation
    mF <- mutationPermutationSwap
    # recombination
    rF <- recombinationPermutationCycleCrossover
    #creation
    cF <- function() sample(n)
    # start optimization
    res <- optimCEGO(x = NULL,
                     fun = fun,
                     control = list(creationFunction=cF,
                                     distanceFunction = dist,
                                     optimizerSettings=list(budget=100,popsize=10,
                                                            mutationFunction=mF,
                                                            recombinationFunction=rF),
                                     evalInit=mi,budget=budget,targetY=0,verbosity=1,
                                     model=modelKriging,
                                     vectorized=FALSE))
    print(res)
    return(list(res$xbest, res$ybest, do.call(rbind, res$x), res$y))
    }
    """

    # with open('myfunc.r', 'r') as f:
    #     rstring = f.read()
    rcode = STAP(rstring, "rcode")
    lop_r_fitness = lop.make_r_fitness()
    best_x, best_fitness, x , y = rcode.my_cego(lop_r_fitness,
                                                dist = r_kendallTau,
                                                n = lop.n,
                                                mi = mi,
                                                budget = budget)
    x = np.asarray(x)
    y = np.asarray(y)
    print(x)
    print(y)
    best_x = np.asarray(best_x)
    best_fitness = np.asarray(best_fitness)[0]
    print(f'best: {best_x}\nbest_fitness: {best_fitness}')


np.random.seed(42)
lop = LOP(10,100, phi=0.9)
y = runCEGO(lop, mi = 10, budget = 15)
