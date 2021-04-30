import numpy as np
import itertools as it
from mallows_model import phi2theta, theta2phi
from scipy.optimize import linear_sum_assignment

def dist_at_uniform(n): return n

def distance(a, b):
    return len(a) - np.sum(a==b)

def weighted_median(sample, ws):
    return uHungarian(sample, ws)

def uHungarian(sample, ws):
    m,n   = sample.shape
    wmarg = np.zeros((n,n))
    for i in range(n): # TODO incremental
      for j in range(n):
        freqs = (sample[:,i]==j)
        wmarg[i,j] = (freqs * ws).sum()
    row_ind, col_ind  = linear_sum_assignment( -wmarg )
    # import matplotlib.pyplot as plt
    # import seaborn as sns
    # if len(ws)%40==0 or len(ws)==10:
    #     print(np.around(wmarg,2))
    #     sns.heatmap(-wmarg)
    #     plt.show()
    #     print(col_ind)
    #     print(np.around(ws,2))
    return col_ind


def sample(m,n, phi,sigma0): # INTERFACE
    sample = np.zeros((m,n))
    theta = phi2theta(phi)

    facts_ = np.array([1,1]+[0]*(n-1),dtype=np.float) # TODO precompute
    deran_num_ = np.array([1,0]+[0]*(n-1),dtype=np.float)
    for i in range(2,n+1):
        facts_[i] = facts_[i-1] * i
        deran_num_[i] = deran_num_[i-1]*(i-1) + deran_num_[i-2]*(i-1);
    hamm_count_ = np.array([ deran_num_[d]*facts_[n] / (facts_[d] * facts_[n - d]) for d in range(n+1)],dtype=np.float)
    probsd = np.array([hamm_count_[d] * np.exp(-theta * d) for d in range(n+1)],dtype=np.float)

    for m_ in range(m):
        target_distance = np.random.choice(n+1,p=probsd/probsd.sum())
        sample[m_,:] = sample_at_dist(n,target_distance, sigma0)
    # print("theta",  round(theta,3), round(expected_dist(n,phi)),target_distance)
    return sample


def sample_at_dist(n,unfixed_points_num, sigma0=None):
    # generate a permutation with fixed_points_num (at hamming distance n-fixed_points_num)
    if sigma0 is None: sigma0=np.arange(n)
    sigma = np.zeros(n)-1
    fixed_points = np.random.choice(n,n-unfixed_points_num, replace=False)
    sigma[fixed_points] = fixed_points
    unfix = np.setdiff1d(np.arange(n),fixed_points)
    unfix = np.random.permutation(unfix)
    for i in range(len(unfix)-1):
      sigma[unfix[i]] = unfix[i+1]
    if len(unfix) > 0 : sigma[unfix[-1]] = unfix[0]
    return sigma[sigma0]



def expected_dist(n,phi):
    facts_ = np.array([1,1]+[0]*(n-1),dtype=np.float) # TODO precompute
    # deran_num_ = np.array([1,0]+[0]*(n-1),dtype=np.float)
    for i in range(2,n+1):
        facts_[i] = facts_[i-1] * i
    x_n_1 , x_n= 0,0
    theta = phi2theta(phi)
    for k in range(n+1):
        aux = (np.exp(theta)-1)**k / facts_[k]
        x_n += aux
        if k<n : x_n_1 += aux
    return (n * x_n - x_n_1 * np.exp( theta )) / x_n

#     double Hamming::expectation(double theta){
#     double x_n = 0 , x_n_1 = 0, aux = 0 ;
#     for (int k = 0 ; k <= n_ ; k ++){
#         aux = pow (exp(theta )-1, k) / facts_[ k ];
#         x_n += aux;
#         if (k < n_ )
#             x_n_1 += aux ;//pow (exp(theta )-1, k) / facts_[ k ];
#     }
#     return (double )(n_ * x_n - x_n_1 * exp( theta )) / x_n;
# }



# FIXME: COPY PASTEQ!!!! The comment says it searches for theta, but the function is called find_phi
# These two functions search for theta for different E[D].
# from 0 < E[D] < 1 (large theta) to E_0[D] (theta=0)
# copy paste in hamming
def find_phi(n, dmin, dmax): #NO
    imin, imax = 0.0, 1.0
    iterat = 0
    while iterat < 500:
        med = (imax + imin) / 2
        # FIXME: Here we convert phi2theta, but expected_dist_MM then convert theta to phi???
        d = expected_dist(n, med)
        #print(imin, imax, med, d,imin==imax)
        if d < dmin : imin = med
        elif d > dmax: imax = med
        else: return med
        iterat  += 1
    # FIXME: Is there a default?
    assert False
