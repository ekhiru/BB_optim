from imp import reload
import numpy as np
import mallows_kendall as mk
from scipy.spatial import distance
import pandas as pd


def binary_search_rho(w, ratio_samples_learn, weight_mass_learn,
                      # 0 <= w_i <= 1, w is sorted increasingly,
                      rho_ini=1, rho_end=0, tol=0.001):
  w = np.asarray(w)
  assert np.all(w >= 0.0)
  assert np.all(w <= 1.0)

  # If pos is None we take the largest 4th.
  # Find the rho s.t. the largest 25%(ratio_samples) of the weights  (rho**ws) take the 0.9(weight_mass) of the total ws.  rho^w[:pos] = 0.9*rho^w
  # codes as a recursive binary search in (0,1)
  pos = int(len(w) * ratio_samples_learn)
  rho_med = rho_ini + (rho_end - rho_ini) / 2
  try:
      acum = np.cumsum(rho_med ** w)
      a = acum[pos]
      b = acum[-1]
      if b < tol: # If b is very small, all values are equal, the value of rho does not matter. Let's return 1.0
        return 1.0
      if abs(a / b - weight_mass_learn) < tol:
          return rho_med
      if a / b > weight_mass_learn:
          mid, last = rho_ini, rho_med
      else:
          mid, last = rho_med, rho_end
      return binary_search_rho(w, ratio_samples_learn, weight_mass_learn, mid, last)
  except: # MANUEL: How can the above fail?
       print(w)
       pos = int(len(w) * ratio_samples_learn)
       print(pos,len(w),ratio_samples_learn)
       rho_med = rho_ini + (rho_end - rho_ini) / 2
       acum = np.cumsum(rho_med ** w)
       a = acum[pos]
       b = acum[-1]
       print("Error in binary_search_rho:",
             a,b,a/b,weight_mass_learn,pos,w, rho_med,ratio_samples_learn, weight_mass_learn,rho_ini,rho_end)
       raise

def get_expected_distance(iterat, n, budget):
  # MANUEL: Should this be Kendall max dist?
  N = (n - 1) * n / 2
  f_ini, f_end = N / 4, 1
  iter_decrease = budget - 10 # MANUEL: Why 10?
  jump = (f_ini - f_end) / iter_decrease
  a = f_ini - jump * iterat
  return max(a, f_end)


def UMM(instance, seed, budget,
        m_ini, budgetMM,
        ratio_samples_learn,
        weight_mass_learn, eval_ranks):

    np.random.seed(seed)

    if eval_ranks: # If True, the objective function works with ranks
      # FIXME: Do we really need this lambda?
      f_eval = lambda p: instance.fitness(p)
    else: # Otherwise, it works with orders
      f_eval = lambda p: instance.fitness(np.argsort(p))

    n = instance.n
    sample = [ np.random.permutation(np.arange(n)) for _ in range(m_ini)]
    fitnesses = [f_eval(perm) for perm in sample]
    # ['rho','phi_estim','phi_sample','Distance']
    res = [ [np.nan, np.nan, np.nan,
             instance.distance_to_best(perm)] for perm in sample]

    for m in range(budget - m_ini):
        ws = np.asarray(fitnesses).copy()
        # FIXME: For maximization, this need to be changed.
        ws = ws - ws.min()
        # FIXME: Handle if ws.max() == 0.
        ws = ws / ws.max()
        co = ws.copy()
        co.sort()
        rho = binary_search_rho(co, ratio_samples_learn, weight_mass_learn)
        # print(fitnesses)
        # print(ws)
        ws = rho ** ws #MINIMIZE
        # print(ws)
        # ws = rho ** (1-ws) #MAXIMIZE
        # print(ws,co[:int(len(co)/4)].sum(),co.sum())

        borda = mk.uborda(np.array(sample), ws)
        phi_estim = mk.u_phi(sample, borda, ws)
        expected_dist = get_expected_distance(m, n, budget)
        phi_sample = mk.find_phi(n, expected_dist, expected_dist + 1)
        perms = mk.samplingMM(budgetMM, n, phi=phi_sample, k=None)
        #perm = perm[borda]
        # Transforms from sampling space to Borda space.
        perms = [perm[borda] for perm in perms]
        dists = distance.cdist(perms, sample, metric=mk.kendallTau)
        # MANUEL: We probably do not need to sort, just find the min per axis=1.
        dists = np.sort(dists, axis=1)
        indi = np.argmax(dists[:, 0]) #index of the perm with the farthest closest permutation. Maximizes the min dist to the sample
        perm = perms[indi]
        # FIXME: This should already be an array of int type.
        perm = np.asarray(perm, dtype='int')
        sample.append(perm)
        fitnesses.append(f_eval(perm))
        # print(f"UMM: eval={m}\tF={fitnesses[-1]}\tbest_known={instance.best_fitness}")
        # print(fitnesses,ws)

        # This is only used for reporting stats.
        res.append([rho, phi_estim, phi_sample, instance.distance_to_best(borda)])
    df = pd.DataFrame(res, columns=['rho','phi_estim','phi_sample','Distance'])
    df['Fitness'] = fitnesses
    df['x'] = sample
    df['m_ini'] = m_ini
    df['seed'] = seed
    df['budget'] = budget
    df['budgetMM'] = budgetMM
    df['ratio_samples_learn'] = ratio_samples_learn
    df['weight_mass_learn'] = weight_mass_learn
    df['eval_ranks'] = eval_ranks
    return df
