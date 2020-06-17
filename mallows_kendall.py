import numpy as np
import itertools as it
import scipy as sp
#fit MM

# estas dos funciones valen para buscar los thetas para diferentes E[d].
# desde 0<E[d]<1 (theta alto) hasta E_0[d] (theta 0)
def find_phi_n(n, bins): #NO
    ed, phi_ed = [], []
    #expected dist en la uniforma para este n
    ed_uniform = (n*(n-1)/2)/2 # dist max div 2. AKA
    # ed_uniform = np.mean([expected_dist_MM(n,-0.01), expected_dist_MM(n,0.01)])
    # if k is not None : ed_uniform = np.mean([expected_V(n,theta=None, phi=0.99999,k=k).sum(), expected_V(n,theta=None, phi=1.00000001,k=k).sum()])
    # print(ed_uniform, "find")
    for dmin in np.linspace(0,ed_uniform-1,bins):
        ed.append(dmin)
        phi_ed.append(find_phi(n, dmin, dmin+1))
    return ed, phi_ed
def find_phi(n, dmin, dmax): #NO
    imin, imax = np.float64(0),np.float64(1)
    iterat = 0
    while iterat < 500:
        med = imin + (imax-imin)/2
        d = expected_dist_MM(n,phi2theta(med))#mk.phi2theta(med)
        #print(imin, imax, med, d,imin==imax)
        if d < dmax and d > dmin: return med
        elif d < dmin : imin = med
        elif d > dmax : imax = med
        iterat  += 1

def max_dist(n):
    return int(n*(n-1)/2)

def compose(s,p):
    return np.array(s[p])
def compose_partial(partial,full):#s is partial
    return [partial[i] if not np.isnan(i) else np.nan for i in full]
def inverse_partial(sigma):
    inv = np.array([np.nan]*len(sigma))
    for i,j in enumerate(sigma):
        if not np.isnan(j):
            inv[int(j)] = i
    return inv
def inverse(s):
    return np.argsort(s)

def alpha2beta(alpha,k):
    inv = np.array([np.nan]*len(alpha))
    for i,j in enumerate(alpha[:k]): inv[int(j)] = i
    return inv


def borda(rankings):
    consensus =  np.argsort(np.argsort(rankings.sum(axis=0))) #borda
    return consensus

def borda_partial(rankings):
    borda = np.argsort(np.argsort(np.nanmean(rankings, axis=0))).astype(float)
    mask = np.isnan(rankings).all(axis=0)
    borda[mask]=np.nan
    return borda


def check_theta_phi(theta, phi):
    if not ((phi is None) ^ (theta is None)): print("KAKA, pon valores")
    if phi is None and type(theta)!=list: phi = theta2phi(theta)
    if theta is None and type(phi)!=list: theta = phi2theta(phi)
    if phi is None and type(theta)==list: phi = [theta2phi(t) for t in theta]
    if theta is None and type(phi)==list: theta = [phi2theta(p) for p in phi]
    return theta, phi

def expected_dist_MM(n,theta=None, phi=None):
    theta, phi = check_theta_phi(theta, phi)
    expected_dist = n * np.exp(-theta) / (1-np.exp(-theta)) - np.sum([j * np.exp(-j*theta) / (1 - np.exp(-j*theta))  for j in range(1,n+1)])
    return expected_dist
def variance_dist_MM(n,theta=None, phi=None):
    theta, phi = check_theta_phi(theta, phi)
    variance = (phi*n)/(1-phi)**2 - np.sum([(pow(phi,i) * i**2)/(1-pow(phi,i))**2  for i in range(1,n+1)])
    return variance
def expected_V(n,theta=None, phi=None,k=None):#txapu integrar
    theta, phi = check_theta_phi(theta, phi)
    if k is None: k = n-1
    if type(theta)!=list: theta = [theta]*k
    expected_v = [np.exp(-theta[j]) / (1-np.exp(-theta[j])) - (n-j) * np.exp(-(n-j)*theta[j]) / (1 - np.exp(-(n-j)*theta[j]))  for j in range(k)]
    return np.array(expected_v)
def variance_V(n,theta=None, phi=None,k=None):#txapu integrar es posibe q solo fuciones con MM
    theta, phi = check_theta_phi(theta, phi)
    if k is None: k = n-1
    if type(phi)!=list: phi = [phi]*k
    var_v = [phi[j]/(1-phi[j])**2 - (n-j)**2 * phi[j]**(n-j) / (1-phi[j]**(n-j))**2 for j in range(k)]
    return np.array(var_v)

def psiMM(n,theta=None,phi=None):
    if theta is not None: return np.prod([(1-np.exp(-theta*j))/(1-np.exp(-theta)) for j in range(2,n+1)])
    if phi is not None:  return np.prod([(1-np.power(phi,j))/(1-phi) for j in range(2,n+1)])
    theta, phi = check_theta_phi(theta, phi) #por si acaso
    #np.array([(1 - np.exp(( - n + i )*(theta)))/(1 - np.exp( -theta)) for i in range(n-1)])

def prob_mode(n, theta):
    #theta as array
    psi = np.array([(1 - np.exp(( - n + i )*(theta[i])))/(1 - np.exp( -theta[i])) for i in range(n-1)])
    return np.prod(1.0/psi)

def prob(n, theta, dist):
    psi = np.array([(1 - np.exp(( - n + i )*(theta)))/(1 - np.exp( -theta)) for i in range(n-1)])
    psi = np.prod(psi)
    return np.exp(-theta*dist) / psi

def prob_sample(perms,sigma,theta=None,phi=None):
    m,n = perms.shape
    theta, phi = check_theta_phi(theta, phi)
    psi = np.array([(1 - np.exp(( - n + i )*(theta)))/(1 - np.exp( -theta)) for i in range(n-1)])
    psi = np.prod(psi)
    return np.array([np.exp(-theta*kendallTau(perm, sigma)) / psi for perm in perms])

def fit_MM(rankings, s0=None): #returns sigma, phi
    m , n = rankings.shape
    if s0 is None: s0 = np.argsort(np.argsort(rankings.sum(axis=0))) #borda
    dist_avg = np.mean(np.array([kendallTau(s0, perm) for perm in rankings]))
    try:
        theta = sp.optimize.newton(mle_theta_mm_f, 0.01, fprime=mle_theta_mm_fdev, args=(n, dist_avg), tol=1.48e-08, maxiter=500, fprime2=None)
    except:
        if dist_avg == 0.0: return s0, np.exp(-5)#=phi
        print("error. fit_mm. dist_avg=",dist_avg, dist_avg == 0.0)
        print(rankings)
        print(s0)
        raise
    # theta = - np.log(phi)
    return s0, np.exp(-theta)#=phi

def fit_MM_phi(n, dist_avg): #returns sigma, phi
    try:
        theta = sp.optimize.newton(mle_theta_mm_f, 0.01, fprime=mle_theta_mm_fdev, args=(n, dist_avg), tol=1.48e-08, maxiter=500, fprime2=None)
    except:
        if dist_avg == 0.0: return s0, np.exp(-5)#=phi
        print("error. fit_mm. dist_avg=",dist_avg, dist_avg == 0.0)
        print(rankings)
        print(s0)
        raise
    # theta = - np.log(phi)
    return np.exp(-theta)

def theta2phi(theta):
    return np.exp(-theta)
def phi2theta(phi):
    return - np.log(phi)

def mle_theta_mm_f(theta, n, dist_avg):
    aux = 0
    for j in range(1,n):
        k = n - j + 1
        aux += (k * np.exp(-theta * k))/(1 - np.exp(-theta * k))
    aux2 = (n-1) / (np.exp( theta ) - 1) - dist_avg
    return aux2 - aux

def mle_theta_mm_fdev(theta, n, dist_avg):
    aux = 0
    for j in range(1,n):
        k = n - j + 1
        aux += (k * k * np.exp( -theta * k ))/pow((1 - np.exp(-theta * k)) , 2 )
    aux2 = (- n + 1) * np.exp( theta ) / pow ((np.exp( theta ) - 1) , 2 )
    # print(theta)
    return aux2 + aux

def likelihood_mm(perms, s0, theta):
    m,n = perms.shape
    psi = 1.0 / np.prod([(1-np.exp(-theta*j))/(1-np.exp(-theta)) for j in range(2,n+1)])
    probs = np.array([np.log(np.exp(-kendallTau(s0, perm)*theta)/psi) for perm in perms])
    # print(probs,m,n)
    return probs.sum()

def samplingMM(m,n,theta=None, phi=None, k=None):
    # k return partial orderings
    theta, phi = check_theta_phi(theta, phi)
    if k==n:k=None
    return samplingGMM(m,[theta]*(n-1),topk=k)

def samplingGMM(m,theta, topk=None):
    #  returns RANKINGS!!!!!!!*****
    n = len(theta)+1
    if topk is None or topk == n: k = n-1
    else: k = topk
    psi = [(1 - np.exp(( - n + i )*(theta[ i ])))/(1 - np.exp( -theta[i])) for i in range(k)]
    vprobs = np.zeros((n,n))
    for j in range(k): #range(n-1):
        vprobs[j][0] = 1.0/psi[j]
        for r in range(1,n-j):
            vprobs[j][r] = np.exp( -theta[j] * r ) / psi[j]#vprobs[j][ r - 1 ] + np.exp( -theta[j] * r ) / psi[j]
    sample = []
    vs = []
    for samp in range(m):
        v = [np.random.choice(n,p=vprobs[i,:]) for i in range(k)] # v = [np.random.choice(n,p=vprobs[i,:]/np.sum(vprobs[i,:])) for i in range(k)]
        #vs.append(v)
        #print(v, np.sum(v))
        # print(v, topk)
        if topk is None: v += [0] # la fun discordancesToPermut necesita, len(v)==n
        ranking = v2ranking(v, n)#discordancesToPermut(v,list(range(n)))
        # if topk is not None :
        #     ranking = np.concatenate([ranking, np.array([np.nan]*(n-topk))])
        sample.append(ranking)
    return sample

def ranking2v(perm):
    n = len(perm)
    return np.array([np.sum([perm[i]<perm[j] for i in range(j+1,n)], dtype=int) for j in range(n)])

def ranking2vinv(perm):
    inv = np.argsort(perm)
    n = len(perm)
    return np.array([np.sum([inv[i]<inv[j] for i in range(j+1,n)], dtype=int) for j in range(n)])

def v2ranking(v, n): ##len(v)==n, last item must be 0
    # n = len(v)
    rem = list(range(n))
    rank = np.array([np.nan]*n)# np.zeros(n,dtype=np.int)
    # print(v,rem,rank)
    for i in range(len(v)):
        # print(i,v[i], rem)
        rank[i] = rem[v[i]]
        rem.pop(v[i])
    return rank#[i+1 for i in permut];


def discordancesToPermut(indCode, refer):
    print("warning. discordancesToPermut is deprecated. Use function v2ranking")
    return v2ranking(indCode)
    # returns rNKING
    # n = len(indCode)
    # rem = refer[:] #[i for i in refer]
    # ordering = np.zeros(n,dtype=np.int)
    # for i in range(n):
    #     ordering[i] = rem[indCode[i]]
    #     rem.pop(indCode[i])
    # return ordering#[i+1 for i in permut];

def kendallTau(A, B=None):
    # if any partial is B
    if B is None : B = list(range(len(A)))
    n = len(A)
    pairs = it.combinations(range(n), 2)
    distance = 0
    # print("IIIIMNNMNNN",list(pairs),len(A))
    for x, y in pairs:
        #if not A[x]!=A[x] and not A[y]!=A[y]:#OJO no se check B
        a = A[x] - A[y]
        try:
            b = B[x] - B[y]# if discordant (different signs)
        except:
            print("ERROR kendallTau, check b",A, B, x, y)
        # print(b,a,b,A, B, x, y,a * b < 0)
        if (a * b < 0):
            distance += 1
    return distance


def partial_ord2partial_rank(pord,n,k,type="beta"):#NO
    if type=="gamma": val = -1
    if type=="beta": val = k

    # pord is a collection of partial orderings, each of which (1) has len n (2) np.nans for the unspecified places (3) is np.array
    #input partial ordering of the first k items. The first k positions have vals [0,n-1]
    #output partial ranking of the first k ranks. There are k positions have vals [0,k-1]. The rest have val=k (so the kendall dist can be compared)
    prank = []
    # n = len(pord[0])
    # for perm in pord:
    res = np.array([val]*n)
    for i,j in enumerate(pord[~np.isnan(pord)]):
        res[int(j)]=i
    # prank.append(res)
    return np.array(res)

# m'/M segun wolfram -((j - n) e^(j x))/(e^(n x) - e^(j x)) - (j e^x - j - n e^x + n + e^x)/(e^x - 1)
#
def Ewolfram(n,j,x):#NO
    return (-((j - n) * np.exp(j * x))/(np.exp(n* x) - np.exp(j *x)) - (j* np.exp(x) - j - n *np.exp(x) + n + np.exp(x))/(np.exp(x) - 1))
#(E^(x + 2 j x) + E^(x + 2 n x) - E^((j + n) x) (j - n)^2 - E^((2 + j + n) x) (j - n)^2 + 2 E^((1 + j + n) x) (-1 + j^2 - 2 j n + n^2))/((-1 + E^x)^2 (-E^(j x) + E^(n x))^2)
def Vwolfram(n,j,x):#NO
    numer = (np.exp(x + 2* j *x) + np.exp(x + 2* n *x) - np.exp((j + n) *x)*(j - n)**2 - np.exp((2+j+n)* x)* (j - n)**2 + 2 *np.exp((1 + j + n) *x)*(-1 + j**2 - 2* j *n + n**2))
    denom = ((-1 + np.exp(x))**2 *(-np.exp(j *x) + np.exp(n* x))**2)
    return numer/denom


## number of perms at each dist
def num_perms_at_dist(n):
    sk = np.zeros((n+1,int(n*(n-1)/2+1)))
    for i in range(n+1):
        sk[i,0] = 1
    for i in range(1,1+n):
        for j in range(1,int(i*(i-1)/2+1)):
            if j - i >= 0 :
                sk[i,j] = sk[i,j-1]+ sk[i-1,j] - sk[i-1,j-i]
            else:
                sk[i,j] = sk[i,j-1]+ sk[i-1,j]
    return sk.astype(np.uint64)
## random permutations at distance
def random_perm_at_dist(n, dist, sk):
    # param sk is the results of the function num_perms_at_dist(n)
    i = 0
    probs = np.zeros(n+1)
    v = np.zeros(n,dtype=int)
    while i<n and dist > 0 :
        rest_max_dist = (n - i - 1 ) * ( n - i - 2 ) / 2
        if rest_max_dist  >= dist:
            probs[0] = sk[n-i-1,dist]
        else:
            probs[0] = 0
        mi = min(dist + 1 , n - i )
        for j in range(1,mi):
            if rest_max_dist + j >= dist: probs[j] = sk[n-i-1, dist-j]
            else: probs[ j ] = 0
        v[i] = np.random.choice(mi,1,p=probs[:mi]/probs[:mi].sum())
        dist -= v[i]
        i += 1
    return v2ranking(v)











# end
