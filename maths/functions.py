import numpy as np


def ctFet(eta):
    '''Gives cosine of scattering angle given pseudorapidity .

    Parameters
    -----------
    eta : float
        Pseudorapidity.

    Returns
    --------
    float
        Cosine of scattering angle
    '''
    return np.cos(2*np.arctan(np.exp(-eta)))


def thFet(eta):
    '''Gives the scattering angle given pseudorapidity .

    Parameters
    -----------
    eta : float
        Pseudorapidity.

    Returns
    --------
    flooat
        The scattering angle
    '''
    return 2*np.arctan(np.exp(eta))


def etFct(ct):
    '''Gives the pseudorapidity given the cosine of the scattering angle .

    Parameters
    ----------
    ctheta : float
        Cosine of scatterting angle

    Returns
    --------
    float
        pseudorapidity
    '''
    return -np.log(np.tan(np.arccos(ct)/2.))


def lamK(x, y, z):
    '''Returns the Kallen function

    Parameters
    ----------
    x, y, z : floats
        The variables entering the Kallen function

    Returns
    --------
    float
        The Kallen Function
    '''
    return x**2 + y**2 + z**2 - 2*(x*y + y*z + z*x)
