import scipy as sp
import time
import sys
import os
import numpy as np
import pandas as pd
import mallows_kendall as mk
import lop
from lop import LOP
import cego
from imp import reload
reload(cego)
reload(lop)

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


#def get_expected_distance(iter):    #return ini_dist - ini_dist * iter_ratio +2
def get_expected_distance(iterat,n,m):
    N = (n-1)*n/2
    f_ini, f_end = N/4,1
    iter_decrease = m-10
    salto = (f_ini - f_end ) / iter_decrease
    a = f_ini - salto * iterat
    return max(a,f_end)

#get_expected_distance(0, 4)

def binary_search_rho(w,rho_ini=1,rho_end=0, tol=0.001,ratio_samples=0.25,weight_mass=0.9):#0<w_i<1, w is sorted increasingly,
  #if pos is None we take the largest 4th.
  #find the rho s.t. the largest 25%(ratio_samples) of the weights  (rho*ws) take the 0.9(weight_mass) of the total ws.  rho^w[:pos] = 0.9*rho^w
  # codes as a recursive binary search in (0,1)
  pos = int(len(w)*ratio_samples)
  rho_med = rho_ini + (rho_end-rho_ini)/2
  acum = np.cumsum(rho_med**w)
  a = acum[pos]
  b = acum[-1]
  if abs( a/b - weight_mass) < tol: return rho_med
  if a/b > weight_mass: return binary_search_rho(w,rho_ini,rho_med)
  return binary_search_rho(w,rho_med,rho_end)

def solve_one_umm(instance, ms, rep,  m_ini, best_sol,worst_sol,budgetMM):
    res = []
    n = instance.n
    N = (n-1)*n/2
    sample = [np.random.permutation(range(n)) for _ in range(m_ini)]
    fitnesses = [instance.get_fitness(perm) for perm in sample]
    best_known_fit = instance.get_fitness(best_sol)
    worst_known_fit = instance.get_fitness(worst_sol)
    autorho = []
    for m in range(ms):


        # #if m%10 == 9:print("ws, fitnesses",ws, fitnesses)
        # ws = np.asarray(fitnesses.copy())
        # ws = ws-ws.min()
        # ws = ws/ws.max() # 0<=fitness<=1
        # ws = rho**(ws) #MINIMIZE
        # co = ws.copy()
        # co.sort()
        # co = co[::-1]
        # acum = np.cumsum(co)/co.sum()
        # # print(len(acum[acum<0.9])/len(acum), len(acum[acum<0.9]))
        # autorho.append([len(acum[acum<0.9])/len(acum), len(acum[acum<0.9])])

        ws = np.asarray(fitnesses.copy())
        ws = ws-ws.min()
        ws = ws/ws.max()
        co = ws.copy()
        co.sort()
        rho = binary_search_rho(co)

        # print(rho)
        ws = rho**(ws) #MINIMIZE
        # print(ws,co[:int(len(co)/4)].sum(),co.sum())

        borda = mk.uborda(np.array(sample),ws)

        phi_estim = mk.u_phi(sample,borda, ws)
        expected_dist = get_expected_distance(m,n,ms)
        phi_sample = mk.find_phi(n, expected_dist, expected_dist+1)
        #phi_estim = 1 - (m+1)/(ms)
        perms = mk.samplingMM(budgetMM,n, phi=phi_sample, k=None)
        #perm = perm[borda]
        perms = [perm[borda] for perm in perms]
        dists = sp.spatial.distance. cdist (perms, sample, metric=mk.kendallTau)
        dists = np.sort(dists, axis=1)
        indi = np.argmax(dists[:,0]) #index of the perm with the farthest closest permutation. Maximizes the min dist to the sample
        perm = perms[indi]

        sample.append(perm)
        #print(m,mk.kendallTau(borda,best_sol)/N,mk.kendallTau(perm,best_sol)/N,expected_dist)
        fitnesses.append(instance.get_fitness(perm))
        #print(perm,fitnesses[-1])
        fit_normalizado = (fitnesses[-1]- best_known_fit)/ (worst_known_fit - best_known_fit)
        res.append([instance.problem_name,"uMM ",rep,m,rho,fit_normalizado,phi_estim,phi_sample,mk.kendallTau(borda,best_sol)/(n*(n-1)/2),budgetMM])
    df = pd.DataFrame(res, columns=['Problem','Solver','rep','Sample size','rho','Fitness','phi_estim','phi_sample','Distance','budget'])
    # autorho = pd.DataFrame(autorho,columns=['propor','total'])
    # autorho.total.plot()
    # import matplotlib.pyplot as plt
    # plt.show()
    # autorho.propor.plot()
    return df


def run_and_save(n,rep,phi_instance,budgetGA,budgetMM,SLURM_JOB_ID="Local",m_max=400):
    np.random.seed(rep)
    true_sol = list(range(n))
    m_inst = 200
    m_ini = 10
    # rhos=[1e-3,1e-5,1e-7,1e-9]
    # rhos = [1e-9]
    #rhos = [1e-7]
    instance = LOP.generate_synthetic(n, m_inst, phi_instance)
    problem = instance.problem_name
    start_time = time.time()
    df = pd.DataFrame()
    #df = cego.runCEGO(instance = instance, m_ini = m_ini, m = m_max,rep=rep, best_known_sol=true_sol, worst_known_sol=true_sol[::-1], budgetGA=budgetGA)
    df['run_time'] = time.time() - start_time
    #for rho in rhos:
    start_time = time.time()
    dfuMM = solve_one_umm(instance, m_max, rep, m_ini, true_sol,true_sol[::-1],budgetMM)
    dfuMM['run_time'] = time.time() - start_time
    df = pd.concat([df,dfuMM],sort=False)
    df['best_known'] = instance.get_fitness(true_sol)
    df['worst_known'] = instance.get_fitness(true_sol[::-1])
    df['phi_instance'] = phi_instance
    df['budget'] = budgetGA
    df['n'] = n
    df['m_max'] = m_max
    df.to_pickle('pickles/pick'+str(SLURM_JOB_ID)+'.pkl')


if __name__ == '__main__':
    print(sys.argv)
    params = [float(p) for p in sys.argv[1:]]
    print(params)
    [n,rep,phi_instance,budgetGA,budgetMM,SLURM_JOB_ID,m_max] = params
    print("assigned",[n,rep,phi_instance,budgetGA,budgetMM,SLURM_JOB_ID,m_max] )
    run_and_save(int(n),int(rep),float(phi_instance),int(budgetGA),int(budgetMM),int(SLURM_JOB_ID),int(m_max))
