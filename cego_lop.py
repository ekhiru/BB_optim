import time
import sys
import os
import numpy as np
import pandas as pd
import mallows_kendall as mk
import lop
import cego

def solve_one_umm(problem, instance,ms, rho, rep,  m_ini, true_sol):
    res = []
    n = instance.shape[0]
    sample = [np.random.permutation(range(n)) for _ in range(m_ini)]
    fitnesses = [get_fitness(perm, instance,problem) for perm in sample]
    best_known_fit = lop.get_fitness(true_sol, instance, 'LOP')
    for m in range(ms):
        #if m%10 == 9:print("ws, fitnesses",ws, fitnesses)
        ws = np.array(fitnesses.copy())
        ws = ws-ws.min()
        ws = ws/ws.max()
        ws = rho**(1-ws)
        borda = lop.uborda(np.array(sample),ws)
        phi_estim = lop.u_phi(sample,borda, ws)
        expected_dist = lop.get_expected_distance((m+1)/ms,(n-1)*n/4)#initial distance is the expectred at uniformity
        phi_sample = mk.find_phi(n, expected_dist, expected_dist+1)
        #phi_estim = 1 - (m+1)/(ms)
        perm = mk.samplingMM(1,n, phi=phi_sample, k=None)[0]
        perm = perm[borda]
        sample.append(perm)
        fitnesses.append(lop.get_fitness(perm, instance,problem))
        #print(perm,fitnesses[-1])
        res.append([problem,"uMM, rho= "+str(rho),rep,m,rho,fitnesses[-1]/best_known_fit,phi_estim,phi_sample,mk.kendallTau(borda,true_sol)])
    df = pd.DataFrame(res, columns=['Problem','Solver','rep','Sample size','rho','Fitness','phi_estim','phi_sample','Distance'])
    return df


#np.random.seed(42)
#lop = LOP(10,100, phi=0.9)
#y = runCEGO(lop, mi = 10, budget = 15)

def run_and_save(n,rep,phi_instance,budgetGA,SLURM_JOB_ID="Local",m_max=400):
    #param rho no vale!!! OJO
    true_sol = list(range(n))
    m_inst = 200
    m_ini = 10
    rhos = [0.000001,0.00001,0.0001,0.001,0.01,0.1,0.2,0.3]
    problem = "LOP"
    instance = lop.synthetic_LOP(n,m_inst,phi_instance)
    start_time = time.time()
    #df = pd.DataFrame()
    df = cego.runCEGO(instance = lop.LOP(n, instance), m_ini = m_ini, m = m_max,rep=rep, best_known_sol=true_sol, budgetGA=budgetGA)
    df['run_time'] = time.time() - start_time
    for rho in rhos:
        start_time = time.time()
        dfuMM = solve_one_umm(problem,instance,m_max, rho, rep,  m_ini,true_sol)
        dfuMM['run_time'] = time.time() - start_time
        df = pd.concat([df,dfuMM],sort=False)
    df['best_known'] = lop.get_fitness(true_sol, instance,problem)
    df['worst_known'] = lop.get_fitness(true_sol[::-1], instance,problem)
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
