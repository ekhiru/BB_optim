import time
import sys
import os
import numpy as np
import pandas as pd
import mallows_kendall as mk
from lop import LOP
import cego

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

def solve_one_umm(instance, ms, rho, rep,  m_ini, true_sol):
    res = []
    n = instance.n
    sample = [np.random.permutation(range(n)) for _ in range(m_ini)]
    fitnesses = [instance.get_fitness(perm) for perm in sample]
    best_known_fit = instance.get_fitness(true_sol)
    for m in range(ms):
        #if m%10 == 9:print("ws, fitnesses",ws, fitnesses)
        ws = np.asarray(fitnesses.copy())
        ws = ws-ws.min()
        ws = ws/ws.max()
        ws = rho**(1-ws)
        borda = mk.uborda(np.array(sample),ws)
        phi_estim = mk.u_phi(sample,borda, ws)
        expected_dist = mk.get_expected_distance((m+1)/ms,(n-1)*n/4)#initial distance is the expectred at uniformity
        phi_sample = mk.find_phi(n, expected_dist, expected_dist+1)
        #phi_estim = 1 - (m+1)/(ms)
        perm = mk.samplingMM(1,n, phi=phi_sample, k=None)[0]
        perm = perm[borda]
        sample.append(perm)
        fitnesses.append(instance.get_fitness(perm))
        #print(perm,fitnesses[-1])
        res.append([instance.problem_name,"uMM, rho= "+str(rho),rep,m,rho,fitnesses[-1]/best_known_fit,phi_estim,phi_sample,mk.kendallTau(borda,true_sol)])
    df = pd.DataFrame(res, columns=['Problem','Solver','rep','Sample size','rho','Fitness','phi_estim','phi_sample','Distance'])
    return df


def run_and_save(n,rep,phi_instance,budgetGA,SLURM_JOB_ID="Local",m_max=400):
    np.random.seed(rep)
    #param rho no vale!!! OJO
    true_sol = list(range(n))
    m_inst = 200
    m_ini = 10
    rhos = [0.000001,0.00001,0.0001,0.001,0.01,0.1,0.2,0.3]
    instance = LOP.generate_synthetic(n, m_inst, phi_instance)
    problem = instance.problem_name
    start_time = time.time()
    #df = pd.DataFrame()
    df = cego.runCEGO(instance = instance, m_ini = m_ini, m = m_max,rep=rep, best_known_sol=true_sol, budgetGA=budgetGA)
    df['run_time'] = time.time() - start_time
    for rho in rhos:
        start_time = time.time()
        dfuMM = solve_one_umm(instance,m_max, rho, rep,  m_ini,true_sol)
        dfuMM['run_time'] = time.time() - start_time
        df = pd.concat([df,dfuMM],sort=False)
    df['best_known'] = instance.get_fitness(true_sol)
    df['worst_known'] = instance.get_fitness(true_sol[::-1])
    df['phi_instance'] = phi_instance
    df['budgetGA'] = budgetGA
    df['n'] = n
    df.to_pickle('pickles/pick'+str(SLURM_JOB_ID)+'.pkl')


if __name__ == '__main__':
    print(sys.argv)
    params = [float(p) for p in sys.argv[1:]]
    print(params)
    [n,rep,phi_instance,budgetGA,SLURM_JOB_ID] = params
    print("assigned",[n,rep,phi_instance,budgetGA,SLURM_JOB_ID] )
    run_and_save(int(n),int(rep),float(phi_instance),int(budgetGA),int(SLURM_JOB_ID))
