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

def theta_to_phi(theta):
    """This functions converts theta dispersion parameter into phi
        Parameters
        ----------
        theta: float
            Real dispersion parameter
        Returns
        -------
        float
            phi real dispersion parameter
    """
    return np.exp(-theta)

def phi_to_theta(phi):
    """This functions converts phi dispersion parameter into theta
        Parameters
        ----------
        phi: float
            Real dispersion parameter
        Returns
        -------
        float
            theta real dispersion parameter
    """
    return -np.log(phi)

