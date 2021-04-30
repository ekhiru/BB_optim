def select_mm(dist_name):
    if dist_name == 'hamming':
        import mallows_hamming as mh
        return mh
    elif dist_name == 'kendall':
        import mallows_kendall as mk
        return mk
    else:
      raise

import numpy as np

def theta2phi(theta):
    return np.exp(-theta)

def phi2theta(phi):
    return -np.log(phi)

