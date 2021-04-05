# from imp import reload
import numpy as np
import mallows_kendall as mk
import mallows_hamming as mh
from imp import reload
reload(mh)
from scipy.spatial import distance
import pandas as pd


def binary_search_rho(w, ratio_samples_learn, weight_mass_learn, rho_ini=1, rho_end=0, tol=0.001):
  w = np.asarray(w)
   # 0 <= w_i <= 1, w is sorted increasingly,
  assert np.all(w >= 0.0)
  assert np.all(w <= 1.0)

  # If pos is None we take the largest 4th.
  # Find the rho s.t. the largest 25%(ratio_samples) of the weights  (rho**ws) take the 0.9(weight_mass) of the total ws.  rho^w[:pos] = 0.9*rho^w
  # codes as a recursive binary search in (0,1)
  pos = int(len(w) * ratio_samples_learn)
  rho_med = (rho_ini + rho_end) / 2
  # If the interval is very narrow, just return the value.
  if abs(rho_ini - rho_end) < 1E-20:
    return rho_med

  try:
      acum = np.cumsum(rho_med ** w)
      a = acum[pos]
      b = acum[-1]
      # If b is very small, all values are equal, the value of rho does not matter. Let's return 1.0
      if b < tol:
        return 1.0
      # If the differenc eot the target weight_mass is very small, just return.
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
       print(f"binary_search_rho: a={a} b={b} a/b={a/b} wml={weight_mass_learn} rho_med={rho_med} rho_ini={rho_ini} rho_end={rho_end} w={w}")
       raise

def get_expected_distance(iterat, n, budget):
  # MANUEL: Should this be Kendall max dist?
  N = (n - 1) * n / 2
  f_ini, f_end = N / 4, 1
  iter_decrease = budget - 10 # MANUEL: Why 10?
  jump = (f_ini - f_end) / iter_decrease
  a = f_ini - jump * iterat
  return max(a, f_end)

def remove_duplicates(s):
  d = {a.tostring(): a for a in s}
  return list(d.values())

def design_random(m, n):
  """
  m: number of permutations to generate
  n: permutation size"""
  return remove_duplicates([ np.random.permutation(n) for _ in range(m)])

def min_distance(x, s, dist):
  return np.apply_along_axis(dist, -1, np.asarray(s), B=x).min()

def design_maxmindist(m, n, dist = mk.kendallTau, budget = 1000):
  sample = [ np.random.permutation(n) ]
  while len(sample) < m:
    best = np.random.permutation(n)
    best_d = min_distance(best, sample, dist)
    for i in range(budget):
      xnew = np.random.permutation(n)
      xnew_d = min_distance(xnew, sample, dist)
      if xnew_d > best_d:
        best, best_d = xnew, xnew_d
    sample.append(best)
  return remove_duplicates(sample)



def UMM(instance, seed, budget,
        m_ini, budgetMM,
        ratio_samples_learn,
        weight_mass_learn, eval_ranks, init):

    np.random.seed(seed)

    if eval_ranks: # If True, the objective function works with ranks
      # FIXME: Do we really need this lambda?
      f_eval = lambda p: instance.fitness(p)
    else: # Otherwise, it works with orders
      f_eval = lambda p: instance.fitness(np.argsort(p))

    n = instance.n
    if init == "random":
      sample = design_random(m_ini, n)
    elif init == "maxmindist":
      sample = design_maxmindist(m_ini, n)
    else:
      raise f"Invalid init: {init}"

    fitnesses = [f_eval(perm) for perm in sample]
    # ['rho','phi_estim','phi_sample','Distance']
    res = [ [np.nan, np.nan, np.nan,
             instance.distance_to_best(perm)] for perm in sample]

    for m in range(budget - m_ini):
        ws = np.asarray(fitnesses).copy()
        ws = ws - ws.min()
        ws = ws / ws.max()
        co = ws.copy()
        co.sort()
        rho = binary_search_rho(co, ratio_samples_learn, weight_mass_learn)
        ws = rho ** ws #MINIMIZE


        # sigma0 = mk.uborda(np.array(sample), ws) # TODO incremental computing
        sigma0 = mh.uHungarian(np.array(sample), ws) # TODO incremental computing

        phi_estim = mk.u_phi(sample, sigma0, ws)
        expected_dist = get_expected_distance(m, n, budget)
        phi_sample = mk.find_phi(n, expected_dist, expected_dist + 1)

        perms = mk.samplingMM(budgetMM, n, phi=phi_sample, k=None, sigma0=sigma0)
        perm = perms[0]
        # FIXME: This should already be an array of int type.
        perm = np.asarray(perm, dtype='int')
        sample.append(perm)
        fitnesses.append(f_eval(perm))
        # This is only used for reporting stats.
        res.append([rho, phi_estim, phi_sample, instance.distance_to_best(sigma0)])
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
    df['init'] = init
    return df



    # def UMM(instance, seed, budget,
    #         m_ini, budgetMM,
    #         ratio_samples_learn,
    #         weight_mass_learn, eval_ranks, init):
    #
    #     np.random.seed(seed)
    #
    #     if eval_ranks: # If True, the objective function works with ranks
    #       # FIXME: Do we really need this lambda?
    #       f_eval = lambda p: instance.fitness(p)
    #     else: # Otherwise, it works with orders
    #       f_eval = lambda p: instance.fitness(np.argsort(p))
    #
    #     n = instance.n
    #     if init == "random":
    #       sample = design_random(m_ini, n)
    #     elif init == "maxmindist":
    #       sample = design_maxmindist(m_ini, n)
    #     else:
    #       raise f"Invalid init: {init}"
    #
    #     fitnesses = [f_eval(perm) for perm in sample]
    #     # ['rho','phi_estim','phi_sample','Distance']
    #     res = [ [np.nan, np.nan, np.nan,
    #              instance.distance_to_best(perm)] for perm in sample]
    #
    #     for m in range(budget - m_ini):
    #         ws = np.asarray(fitnesses).copy()
    #         # FIXME: We could use rankings for invariance
    #         # FIXME: For maximization, this need to be changed.
    #         ws = ws - ws.min()
    #         # FIXME: Handle if ws.max() == 0.
    #         ws = ws / ws.max()
    #         co = ws.copy()
    #         co.sort()
    #         rho = binary_search_rho(co, ratio_samples_learn, weight_mass_learn)
    #         ws = rho ** ws #MINIMIZE
    #         borda = mk.uborda(np.array(sample), ws) # TODO incremental computing
    #         phi_estim = mk.u_phi(sample, borda, ws)
    #         expected_dist = get_expected_distance(m, n, budget)
    #         phi_sample = mk.find_phi(n, expected_dist, expected_dist + 1)
    #         perms = mk.samplingMM(budgetMM, n, phi=phi_sample, k=None, sigma0=borda)
    #         perm = perms[0]
    #         # FIXME: This should already be an array of int type.
    #         perm = np.asarray(perm, dtype='int')
    #         sample.append(perm)
    #         fitnesses.append(f_eval(perm))
    #         # This is only used for reporting stats.
    #         res.append([rho, phi_estim, phi_sample, instance.distance_to_best(borda)])
    #     df = pd.DataFrame(res, columns=['rho','phi_estim','phi_sample','Distance'])
    #     df['Fitness'] = fitnesses
    #     df['x'] = sample
    #     df['m_ini'] = m_ini
    #     df['seed'] = seed
    #     df['budget'] = budget
    #     df['budgetMM'] = budgetMM
    #     df['ratio_samples_learn'] = ratio_samples_learn
    #     df['weight_mass_learn'] = weight_mass_learn
    #     df['eval_ranks'] = eval_ranks
    #     df['init'] = init
    #     return df
