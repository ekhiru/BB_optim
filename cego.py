import os
os.environ['RPY2_CFFI_MODE'] = "API" # bug in cffi 1.13.0 https://bitbucket.org/rpy2/rpy2/issues/591/runtimeerror-found-a-situation-in-which-we

import numpy as np
from mallows_kendall import kendallTau
import pandas as pd

from rpy2.robjects.packages import importr
from rpy2.robjects import r as R
from rpy2.robjects import numpy2ri
from rpy2.robjects import FloatVector
from rpy2.robjects.packages import STAP
numpy2ri.activate()
import rpy2.rinterface as ri

# funcion de distancia entre permutaciones
@ri.rternalize
def r_kendallTau(A,B):
    return kendallTau(A,B)
    

# rep : replication number.
def runCEGO(instance, m_ini, m, rep, best_known_sol, budgetGA):
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
                                     evalInit=m_ini,budget=budget,targetY=0,verbosity=1,
                                     model=modelKriging,
                                     vectorized=FALSE))
    print(res)
    return(list(res$xbest, res$ybest, do.call(rbind, res$x), res$y))
    }
    """

    # with open('myfunc.r', 'r') as f:
    #     rstring = f.read()
    rcode = STAP(rstring, "rcode")
    r_fitness = instance.make_r_fitness()
    best_x, best_fitness, x , y = rcode.my_cego(r_fitness,
                                                dist = r_kendallTau,
                                                n = instance.n,
                                                m_ini = m_ini,
                                                budget = m,
                                                seed = rep,
                                                budgetGA = budgetGA)

    x = np.asarray(x)
    y = np.asarray(y)
    #best_x = np.asarray(best_x)
    #best_fitness = np.asarray(best_fitness)[0]
    #print(f'best: {best_x}\nbest_fitness: {best_fitness}')
    df = pd.DataFrame()#columns=['problem','rep','m','rho','fitnesses','phi_estim','phi_sample','dist'])
    best_known_fit = instance.get_fitness(best_known_sol)
    df['Fitness'] = instance.evaluations/best_known_fit
    df['Problem'] = 'LOP'
    df['Solver'] = 'CEGO'
    df['Sample size'] = range(m)
    df['rep'] = rep
#    df['budgetGA'] = budgetGA this must be set for all, including uMM, so that we can filter appropriately
    df['Distance'] = [mk.kendallTau(perm,best_known_sol) for perm in instance.solutions]
    return df

