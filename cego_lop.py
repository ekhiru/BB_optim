import time
import sys
from scipy import optimize
import pandas as pd
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
@ri.rternalize
def kendallTau(A, B):
    n = len(A)
    pairs = it.combinations(range(n), 2)
    distance = 0
    for x, y in pairs:
        a = A[x] - A[y]
        try:
            b = B[x] - B[y]# if discordant (different signs)
        except:
            print("ERROR kendallTau, check b",A, B, x, y)
        if (a * b < 0):
            distance += 1
    return distance

def synthetic_LOP(n, m, phi):
  instance = np.zeros((n,n))
  s = np.array(mk.samplingMM(m,n, phi=phi, k=None))
  for i in range(n):
      for j in range(i+1,n):
          instance[i,j] = (s[:,i]< s[:,j]).sum()
          instance[j,i] = m - instance[i,j]
  return instance

def u_phi(sample,s0, ws):
    m , n = np.array(sample).shape
    #if s0 is None: s0 = np.argsort(np.argsort(rankings.sum(axis=0))) #borda
    dist_avg = np.array([mk.kendallTau(perm, s0) for perm in sample]*ws).sum()/ws.sum() #np.mean(np.array([kendallTau(s0, perm) for perm in rankings]))
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
  def __init__(self, n, instance):
    self.n = n
    self.instance = instance
    self.evaluations = []
    self.solutions = []

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
def solve_one_umm(problem, instance,ms, rho, repe,  m_ini, true_sol):
    res = []
    n = instance.shape[0]
    sample = [np.random.permutation(range(n)) for _ in range(m_ini)]
    fitnesses = [get_fitness(perm, instance,problem) for perm in sample]
    best_known_fit = get_fitness(true_sol, instance, 'LOP')
    for m in range(ms):
        #if m%10 == 9:print("ws, fitnesses",ws, fitnesses)
        ws = np.array(fitnesses.copy())
        ws = ws-ws.min()
        ws = ws/ws.max()
        ws = rho**(1-ws)
        borda = uborda(np.array(sample),ws)
        phi_estim = u_phi(sample,borda, ws)
        expected_dist = get_expected_distance((m+1)/ms,(n-1)*n/4)#initial distance is the expectred at uniformity
        phi_sample = mk.find_phi(n, expected_dist, expected_dist+1)
        #phi_estim = 1 - (m+1)/(ms)
        perm = mk.samplingMM(1,n, phi=phi_sample, k=None)[0]
        perm = perm[borda]
        sample.append(perm)
        fitnesses.append(get_fitness(perm, instance,problem))
        #print(perm,fitnesses[-1])
        res.append([problem,"uMM, rho= "+str(rho),repe,m,rho,fitnesses[-1]/best_known_fit,phi_estim,phi_sample,mk.kendallTau(borda,true_sol)])
    df = pd.DataFrame(res, columns=['Problem','Solver','repe','Sample size','rho','Fitness','phi_estim','phi_sample','Distance'])
    return df

def runCEGO(n,instance, m_ini, m, repe, best_known_sol, budgetGA):
    rstring = """
    library(CEGO)
    my_cego <- function(fun, dist, n, m_ini = 5, budget = 15, seed = 0, budgetGA = 100)
    {
    set.seed(seed)
    # mutation
    mF <- mutationPermutationInterchange # changed as paper mutationPermutationSwap
    # recombination
    rF <- recombinationPermutationCycleCrossover
    #creation
    cF <- function() sample(n)
    # start optimization
    res <- optimCEGO(x = NULL,
                     fun = fun,
                     control = list(creationFunction=cF,
                                     distanceFunction = dist,
                                     optimizerSettings=list(budget=budgetGA,popsize=20,
                                                            mutationFunction=mF,
                                                            recombinationFunction=rF),
                                     evalInit=m_ini,budget=budget,targetY=0,verbosity=0,
                                     model=modelKriging,
                                     vectorized=FALSE))
    #print(res)
    return(list(res$xbest, res$ybest, do.call(rbind, res$x), res$y))
    }
    """

    # with open('myfunc.r', 'r') as f:
    #     rstring = f.read()
    lop = LOP(n, instance)
    rcode = STAP(rstring, "rcode")
    lop_r_fitness = lop.make_r_fitness()
    best_x, best_fitness, x , y = rcode.my_cego(lop_r_fitness,
                                                dist = kendallTau,
                                                n = lop.n,
                                                m_ini = m_ini,
                                                budget = m,
                                                seed = repe,
                                                budgetGA = budgetGA)

    x = np.asarray(x)
    y = np.asarray(y)
    #best_x = np.asarray(best_x)
    #best_fitness = np.asarray(best_fitness)[0]
    #print(f'best: {best_x}\nbest_fitness: {best_fitness}')
    df = pd.DataFrame()#columns=['problem','repe','m','rho','fitnesses','phi_estim','phi_sample','dist'])
    best_known_fit = get_fitness(best_known_sol, instance, 'LOP')
    df['Fitness'] = lop.evaluations/best_known_fit
    df['Problem'] = 'LOP'
    df['Solver'] = 'CEGO'
    df['Sample size'] = range(m)
    df['repe'] = repe
#    df['budgetGA'] = budgetGA this must be set for all, including uMM, so that we can filter appropriately
    df['Distance'] = [mk.kendallTau(perm,best_known_sol) for perm in lop.solutions]
    return df


#np.random.seed(42)
#lop = LOP(10,100, phi=0.9)
#y = runCEGO(lop, mi = 10, budget = 15)

def run_and_save(n,repe,phi_instance,budgetGA,SLURM_JOB_ID="Local",m_max=400):
    #param rho no vale!!! OJO
    true_sol = list(range(n))
    m_inst = 200
    m_ini = 10
    rhos = [0.000001,0.00001,0.0001,0.001,0.01,0.1,0.2,0.3]
    problem = "LOP"
    instance = synthetic_LOP(n,m_inst,phi_instance)
    start_time = time.time()
    df = pd.DataFrame()
    df = runCEGO(n,instance, m_ini = m_ini, m = m_max,repe=repe, best_known_sol=true_sol, budgetGA=budgetGA)
    df['run_time'] = time.time() - start_time
    for rho in rhos:
        start_time = time.time()
        dfuMM = solve_one_umm(problem,instance,m_max, rho, repe,  m_ini,true_sol)
        dfuMM['run_time'] = time.time() - start_time
        df = pd.concat([df,dfuMM],sort=False)
    df['best_known'] = get_fitness(true_sol, instance,problem)
    df['worst_known'] = get_fitness(true_sol[::-1], instance,problem)
    df['phi_instance'] = phi_instance
    df['budgetGA'] = budgetGA
    df['n'] = n
    df.to_pickle('pickles/pick'+str(SLURM_JOB_ID)+'.pkl')


if __name__ == '__main__':
    print(sys.argv)
    params = [float(p) for p in sys.argv[1:]]
    print(params)
    [n,repe,phi_instance,budgetGA,SLURM_JOB_ID] = params
    print("assigned",[n,repe,phi_instance,budgetGA,SLURM_JOB_ID] )
    run_and_save(int(n),int(repe),float(phi_instance),int(budgetGA),int(SLURM_JOB_ID))
