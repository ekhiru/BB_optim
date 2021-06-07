import numpy as np
from mallows_model import select_mm
from scipy.spatial import distance
import matplotlib.pyplot as plt
import pandas as pd

def is_duplicated(perm, sample):
  for p in sample:
    if np.array_equal(perm, p):
      return True
  return False

def remove_duplicates(s):
  d = {a.tostring(): a for a in s}
  return list(d.values())

def design_random(m, n):
  """
  m: number of permutations to generate
  n: permutation size"""
  return remove_duplicates([ np.random.permutation(n) for _ in range(m)])

def min_distance(x, s, distance):
  return np.apply_along_axis(distance, -1, np.asarray(s), b=x).min()

def design_maxmindist(m, n, distance, budget = 1000):
  sample = [ np.random.permutation(n) ]
  while len(sample) < m:
    best = np.random.permutation(n)
    best_d = min_distance(best, sample, distance)
    for i in range(budget):
      xnew = np.random.permutation(n)
      xnew_d = min_distance(xnew, sample, distance)
      if xnew_d > best_d:
        best, best_d = xnew, xnew_d
    sample.append(best)
  return remove_duplicates(sample)


def get_expected_distance_at_iteration_t(dist_at_uniform, budget,scalfun='log'):
  # Uniform distance.
  dist_at_uniform -= 2 # to avoid numerical errors

  if scalfun == 'log':
      # MANUEL: Shouldn't this start at zero? so we get e**0 = 1
      # MANUEL: Why not use logspace(np.log(dist_at_uniform), 0) to get the list already reversed?
      dists = np.logspace(.1, np.log(dist_at_uniform),num=budget,base=np.e)[::-1]
  elif scalfun == 'linear':
      # MANUEL: Same here, why not linspace(dist,1) to get the list reversed?
      dists = np.linspace(1,dist_at_uniform,num=budget)[::-1]
  # print("get_expected_distance_at_iteration_t (n)", dists[iterat],n)
  elif scalfun == 'median':
      dists = budget*[0]

  return dists


def get_rho_at_iteration_t (iterat, budget, scalfun='exp'):
  # -1 is uniform, -20 may need some justification.
  # MANUEL: This one is not base e, but the one above is. It seems a bit of inconsistency.
  return np.logspace(-1,-20,budget)[iterat]


def UMM(instance, seed, budget,
        m_ini,
        eval_ranks, init, dist_name, scalfun_learning, scalfun_sampling):

    np.random.seed(seed)

        # print(sample)

    if eval_ranks: # If True, the objective function works with ranks
      # FIXME: Do we really need this lambda?
      f_eval = lambda p: instance.fitness(p)
    else: # Otherwise, it works with orders
      f_eval = lambda p: instance.fitness(np.argsort(p))

    mm = select_mm(dist_name)
    n = instance.n


    import mallows_hamming as mh
    best_manual = np.array([ 10,  7, 26,  3, 19, 13, 22, 24, 29, 15, 27, 11, 17, 21, 18,  8, 12,  5, 16,  0,  9,  2, 20, 14, 25, 23,  6, 28,  1,  4])
    print("meadin onf bests",f_eval(best_manual))
    print("hole1, distance vs fit")
    resplot = []
    for dis in range(2,n+1):
        for _ in range (10):
            resplot.append([f_eval(mh.sample_at_dist(n,dis, sigma0=best_manual)),dis])
    resplot.append([f_eval(np.arange(n)),0])
    print(f_eval(best_manual))
    resplot.append([f_eval(np.arange(n)[::-1]),n])
    print(f_eval(best_manual[::-1]))
    df = pd.DataFrame(resplot,columns=['fit','dis'])
    df.plot.scatter(y='fit', x = 'dis')
    plt.show()


    if init == "random":
      sample = design_random(m_ini, n)
      sample = [np.arange(n) for _ in range(2)] +[ best_manual for _ in range(3)] +[np.random.permutation(n) for _ in range(5)]
      print("OJO,  in UMM.py - manual initielization")
    elif init == "maxmindist":
      sample = design_maxmindist(m_ini, n, distance = mm.distance)
    else:
      raise f"Invalid init: {init}"
    expected_dists = get_expected_distance_at_iteration_t(mm.dist_at_uniform(n), budget, scalfun=scalfun_sampling)


    fitnesses = [f_eval(perm) for perm in sample]
    res = [ [np.nan, np.nan, instance.distance_to_best(perm, mm.distance)] for perm in sample]

    for m in range(budget - m_ini):
        ws = np.asarray(fitnesses).copy()
        ws = ws - ws.min()
        ws = ws / ws.max()
        rho = get_rho_at_iteration_t (m, budget, scalfun_learning)
        ws = rho ** ws # Minimize
        # Get weighted median
        sigma0 = mm.weighted_median(np.array(sample), ws) # TODO incremental computing
        # Update sampling variance decreasing SCHEME
        expected_dist = expected_dists[ m ] # get the expected this at this iteration $m$
        phi_sampling = mm.find_phi(n, expected_dist, expected_dist + 1)
        while True:
          # SAMPLE 1 PERM
          perm = mm.sample(1, n, phi=phi_sampling, s0=sigma0)[0]
          # FIXME: This should already be an array of int type.
          perm = np.asarray(perm, dtype='int')
          # Sample again if the permutation has already been evaluated.
          if not is_duplicated(perm, sample):
            break

        for p in sample:
          assert not np.array_equal(perm, p), f"{perm} found in sample:\n {sample}"

        sample.append(perm)
        fitnesses.append(f_eval(perm))
        theta = -np.log(phi_sampling)
        # print("prob mode , umode", mm.prob_mode(n,theta), mm.proba(range(n),sigma0,theta))
        extra_info = [mm.prob(sigma0,sigma0,theta), mm.prob(np.arange(n),sigma0,theta), mm.distance(perm,sigma0),f_eval(range(n)),f_eval(sigma0)]
        extra_info_cols = ['prob_mode', 'prob_id', 'dist_median_perm','feval_id','feval_median']

        # This is only used for reporting stats.
        # MANUEL: Why sigma0, shouldn't this be perm???
        res.append([rho, phi_sampling, instance.distance_to_best(sigma0, mm.distance),expected_dist]+extra_info)

    df = pd.DataFrame(res, columns=['rho','phi_sampling','Distance','expected_dist']+extra_info_cols)
    df['Fitness'] = fitnesses
    df['x'] = [ ' '.join(map(str,s)) for s in sample ]
    df['m_ini'] = m_ini
    df['seed'] = seed
    df['budget'] = budget
    df['eval_ranks'] = eval_ranks
    df['init'] = init
    df['dist_name'] = dist_name
    df['scalfun_learning'] = scalfun_learning
    df['scalfun_sampling'] = scalfun_sampling
    return df
