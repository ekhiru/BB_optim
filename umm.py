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
  #if pos is None we take the largest 4th.
  #find the rho s.t. the largest 25%(ratio_samples) of the weights  (rho**ws) take the 0.9(weight_mass) of the total ws.  rho^w[:pos] = 0.9*rho^w
  # codes as a recursive binary search in (0,1)
  pos = int(len(w) * ratio_samples_learn)
  # print(pos,len(w),ratio_samples_learn)
  rho_med = rho_ini + (rho_end - rho_ini) / 2
  acum = np.cumsum(rho_med ** w)
  a = acum[pos]
  b = acum[-1]
  if abs(a/b - weight_mass_learn) < tol:
      return rho_med
  if a/b > weight_mass_learn:
      mid = rho_ini
      last = rho_med
  else:
      mid = rho_med
      last = rho_end
  return binary_search_rho(w, ratio_samples_learn, weight_mass_learn, mid, last)

def get_expected_distance(iterat,n,m):
    N = (n - 1) * n / 2
    f_ini, f_end = N / 4, 1
    iter_decrease = m - 10
    salto = (f_ini - f_end ) / iter_decrease
    a = f_ini - salto * iterat
    return max(a, f_end)


def UMM(instance,
        seed, budget,
        m_ini,
        budgetMM,
        ratio_samples_learn,
        weight_mass_learn):
    np.random.seed(seed)
    n = instance.n
    sample = [np.random.permutation(range(n)) for _ in range(m_ini)]
    fitnesses = [instance.fitness(perm) for perm in sample]
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
        #phi_estim = 1 - (m+1)/(budget)
        perms = mk.samplingMM(budgetMM, n, phi=phi_sample, k=None)
        #perm = perm[borda]
        # Transforms from sampling space to board space.
        perms = [perm[borda] for perm in perms]
        dists = distance.cdist(perms, sample, metric=mk.kendallTau)
        # FIXME: We probably do not need to sort, just find the min per axis=1.
        dists = np.sort(dists, axis=1)
        indi = np.argmax(dists[:, 0]) #index of the perm with the farthest closest permutation. Maximizes the min dist to the sample
        perm = perms[indi]
        # FIXME: This should already be an array of int type.
        perm = np.asarray(perm, dtype='int')
        sample.append(perm)
        fitnesses.append(instance.fitness(perm))
        # print("Fitness, best, desde UMM.py", fitnesses[-1], instance.best_fitness)
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
    return df
