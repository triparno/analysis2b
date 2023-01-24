import numpy as np
import cons as cn
import kinematics as kine
import amp as amp


def dfdev(ct, mv, Eph, frame):
    '''Returns the jacobian of the delta function.

    Parameters
    ----------
    ct : float
        Cosine of the scattering angle
    mv : float
        mass of the exotic boson
    Eph : float
        Energy of incoming photon
    frame : str
        Frame of computation, 'lab' or 'com'
    Returns
    -------
    float
        float if 'com'
        tuple of two floats if 'lab'
    '''
    s = cn.mp**2 + 2*Eph*cn.mp
    if frame == "com":
        out = 2*np.sqrt(s)
        return out
    elif frame == "lab":
        ev = kine.EDP(ct, mv, Eph, frame)
        out0 = abs((Eph+cn.mp) - ev[0]*Eph*ct/np.sqrt(ev[0]**2-mv**2))
        out1 = abs((Eph+cn.mp) - ev[1]*Eph*ct/np.sqrt(ev[1]**2-mv**2))
        return (2*out0, 2*out1)


def dSidet(et, mv, Eph, frame):
    '''Differential Distribution of cross-section wrt pseudorapidity
    (coupling strength = 1)

    Parameters
    ----------
    et : float
        pseudorapidity
    mv : float
        mass of the exotic boson
    Eph : float
        Energy of incoming photon
    frame : str
        Frame of computation, 'lab' or 'com'
    Returns
    -------
    float
    '''

    ct = np.cos(2*np.arctan(np.exp(-et)))
    jacobian = 1-ct**2
    s = cn.mp**2 + 2*Eph*cn.mp
    fac = 1/(s-cn.mp**2)/8/np.pi
    etph = 1

    ev = kine.EDP(ct, mv, Eph, frame)
    ampl = amp.AmpSq(ct, mv, Eph, frame)

    def out(a, kin, tht):
        return np.heaviside(tht, 0)*a*kin

    if frame == "com":                      # Mandelstam T
        kin = np.sqrt(ev**2-mv**2)/dfdev(ct, mv, Eph, frame)
        tht = np.sqrt(s) - cn.mp - mv
        res = out(ampl, kin, tht)
    elif frame == "lab":
        tht = Eph + cn.mp
        kin0 = np.sqrt(ev[0]**2-mv**2)/dfdev(ct, mv, Eph, frame)[0]
        kin1 = np.sqrt(ev[1]**2-mv**2)/dfdev(ct, mv, Eph, frame)[1]
        res = out(ampl[0], kin0, tht-ev[0]) + out(ampl[1], kin1, tht-ev[1])

    return etph*fac*res*jacobian*cn.GeVInSq2pb


def dSidct(ct, mv, Eph, frame):
    '''Differential Distribution of cross-section wrt cosine of the
    scattering angle (coupling strength = 1)

    Parameters
    ----------
    ct : float
        Cosine of the scattering angle
    mv : float
        mass of the exotic boson
    Eph : float
        Energy of incoming photon
    frame : str
        Frame of computation, 'lab' or 'com'
    Returns
    -------
    float
    '''

    s = cn.mp**2 + 2*Eph*cn.mp
    fac = 1/(s-cn.mp**2)/8/np.pi
    etph = 1

    ev = kine.EDP(ct, mv, Eph, frame)
    ampl = amp.AmpSq(ct, mv, Eph, frame)

    def out(a, kin, tht):
        return np.heaviside(tht, 0)*a*kin

    if frame == "com":
        kin = np.sqrt(ev**2-mv**2)/dfdev(ct, mv, Eph, frame)
        tht = np.sqrt(s) - cn.mp - mv
        res = out(ampl, kin, tht)
    elif frame == "lab":
        tht = Eph + cn.mp
        kin0 = np.sqrt(ev[0]**2-mv**2)/dfdev(ct, mv, Eph, frame)[0]
        kin1 = np.sqrt(ev[1]**2-mv**2)/dfdev(ct, mv, Eph, frame)[1]
        res = out(ampl[0], kin0, tht-ev[0]) + out(ampl[1], kin1, tht-ev[1])

    return etph*fac*res*cn.GeVInSq2pb*np.heaviside(s-mv-cn.mp, 0)
