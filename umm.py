import matplotlib.pyplot as plt

import numpy as np
import mallows_kendall as mk
import mallows_hamming as mh
from imp import reload
reload(mh)
from scipy.spatial import distance
import pandas as pd


# def remove_duplicates(s):
#   d = {a.tostring(): a for a in s}
#   return list(d.values())

def design_random(m, n):
  """
  m: number of permutations to generate
  n: permutation size"""
  return [ np.random.permutation(n) for _ in range(m)]
  # return remove_duplicates([ np.random.permutation(n) for _ in range(m)])

# def min_distance(x, s, dist):
#   return np.apply_along_axis(dist, -1, np.asarray(s), B=x).min()

# def design_maxmindist(m, n, dist = mk.kendallTau, budget = 1000):
#   sample = [ np.random.permutation(n) ]
#   while len(sample) < m:
#     best = np.random.permutation(n)
#     best_d = min_distance(best, sample, dist)
#     for i in range(budget):
#       xnew = np.random.permutation(n)
#       xnew_d = min_distance(xnew, sample, dist)
#       if xnew_d > best_d:
#         best, best_d = xnew, xnew_d
#     sample.append(best)
#   return remove_duplicates(sample)


def get_expected_distance_at_iteration_t(iterat, n, budget, dist_name='kendall',scalfun='log'):
  if dist_name == 'hamming':
      dist_at_uniform = n
  elif dist_name == 'kendall':
      dist_at_uniform = (n - 1) * n / 4
  else: raise

  # Uniform distance.
  dist_at_uniform -= 2 # to avoid numerical errors

  if scalfun == 'log':
      dists = np.logspace(1,np.log(dist_at_uniform),num=budget,base=np.exp(1))[::-1]
  elif scalfun == 'linear':
      dists = np.linspace(1,dist_at_uniform,num=budget)[::-1]
  # print("get_expected_distance_at_iteration_t (n)", dists[iterat],n)
  return dists[iterat]


def get_rho_at_iteration_t (iterat, budget, scalfun='exp'):
  # -1 is uniform, -20 may need some justification.
  return np.logspace(-1,-20,budget)[iterat]


def UMM(instance, seed, budget,
        m_ini,
        eval_ranks,init, dist_name, scalfun_learning, scalfun_sampling):

    np.random.seed(seed)

    if eval_ranks: # If True, the objective function works with ranks
      # FIXME: Do we really need this lambda?
      f_eval = lambda p: instance.fitness(p)
    else: # Otherwise, it works with orders
      f_eval = lambda p: instance.fitness(np.argsort(p))

    n = instance.n
    if init == "random":
      sample = design_random(m_ini, n)
    #   # sample = [np.arange(n) for _ in range(2)] +[np.arange(n)[::-1] for _ in range(1)] +[np.random.permutation(n) for _ in range(7)]
    #   # print("OJO,  in UMM.py - ini to best", sample)
    # elif init == "maxmindist":
    #   sample = design_maxmindist(m_ini, n)
    # else:
    #   raise f"Invalid init: {init}"

    fitnesses = [f_eval(perm) for perm in sample]
    res = [ [np.nan, np.nan, instance.distance_to_best(perm)] for perm in sample]

    if dist_name == 'hamming':
        mmdist = mh
    else: mmdist = mk

    for m in range(budget - m_ini):
        ws = np.asarray(fitnesses).copy()
        ws = ws - ws.min()
        ws = ws / ws.max()
        rho = get_rho_at_iteration_t (m, budget, scalfun_learning)
        ws = rho ** ws # Minimize
        # Get weighted median
        sigma0 = mmdist.weighted_median(np.array(sample), ws) # TODO incremental computing
        # Update sampling variance decreasing SCHEME
        expected_dist = get_expected_distance_at_iteration_t(m, n, budget, dist_name=dist_name, scalfun=scalfun_sampling)
        phi_sample = mmdist.find_phi(n, expected_dist, expected_dist + 1)
        # SAMPLE 1 PERM
        perm = mmdist.sample(1, n, phi=phi_sample, sigma0=sigma0)[0]
        # FIXME: This should already be an array of int type.
        perm = np.asarray(perm, dtype='int')
        sample.append(perm)
        fitnesses.append(f_eval(perm))

        # This is only used for reporting stats.
        res.append([rho, phi_sample, instance.distance_to_best(sigma0)])

    df = pd.DataFrame(res, columns=['rho','phi_sample','Distance'])
    df['Fitness'] = fitnesses
    df['x'] = sample
    df['m_ini'] = m_ini
    df['seed'] = seed
    df['budget'] = budget
    df['eval_ranks'] = eval_ranks
    df['init'] = init
    return df
