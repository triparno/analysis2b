import cons as cn
import numpy as np


def gam_fr(Eph):
    '''Gives the relativistic gamma factor between the com and the lab
    frames.

    Parameters
    ----------
    Eph : float
        Energy of the incoming photon
    Returns
    -------
    gam_fr
        float, the relativistic gamma factor
    '''
    return (Eph+cn.mp)/np.sqrt(cn.mp**2+2*cn.mp*Eph)


def bet_fr(Eph):
    '''Gives the boost between the com and the lab frames.

    Parameters
    ----------
    Eph : float
        Energy of the incoming photon
    Returns
    -------
    bet_fr
        float, the boost
    '''
    return Eph/np.sqrt(cn.mp**2+2*cn.mp*Eph)


def EDP(ct, mv, Eph, frame):
    '''Returns the energy of the DP given the cosine of the scattering
    angle, the mass of the exotic vector, the incoming photon beam
    energy, and the frame of computations.

    Parameters
    ----------
    ct: float
        Cosine of Azimuthal Angle
    mv: float
        Mass of the exotic vector
    Eph: float
        Energy of incoming photon
    frame: str, optional
        Frame of computation, 'lab' or 'com' (default = 'com')

    Returns
    -------
    EDP
        float for com, energy of the outgoing exotic vector
        corresponding to the single root
        OR
        tuple of floats for lac, energy of the outgoing exotic vector
        corresponding to the two roots.
    '''

    s = cn.mp**2 + 2*Eph*cn.mp
    if frame == "com":
        out = (s + mv**2 - cn.mp**2)/np.sqrt(s)/2
        return out
    elif frame == "lab":
        denom = (Eph+cn.mp)**2 - Eph**2*ct**2
        numot = (Eph + cn.mp)*(mv**2 + 2*Eph*cn.mp)
        numin = (
                 (2*Eph*mv*ct)**2
                 + (mv*(mv-2*cn.mp)+2*Eph*(cn.mp-mv))
                 * (mv*(mv+2*cn.mp)+2*Eph*(cn.mp+mv))
                 )
        out1 = (numot - ct*Eph*np.sqrt(numin)) / 2 / denom
        out2 = (numot + ct*Eph*np.sqrt(numin)) / 2 / denom
        return (out1, out2)
