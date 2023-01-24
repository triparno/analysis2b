import cons as cn
import kinematics as kin
import numpy as np


def AmpSq(ct, mv, Eph, frame):
    '''Returns the amplitude---for coupling=1---given the DP mass
    and the cosine of the azimutahl angle.

    Parameters
    ------
    ct: float
        Cosine of the Azimuthal Angle
    mv: float
        Mass of the exotic vector
    Eph: float
        Energy of incoming photon
    frame: str
        Frame of computation
    Returns
    -------
    float
        spin and helicity summed amplitude-squared
    '''

    s = cn.mp**2 + 2*Eph*cn.mp
    faca = cn.aEM * np.pi / 2
    ev = kin.EDP(ct, mv, Eph, frame)

    def trc(t):
        denom = 32/mv**2/(cn.mp**2-s)**2/(cn.mp**2-t)**2
        num0 = (2*(cn.mp**2 - s)*(cn.mp**2 - t)**2
                * (cn.mp**4 + cn.mp**2*s - s*(s + t)))
        num2 = (4*cn.mp**8 + 4*cn.mp**6*(6*s - t)
                - cn.mp**4*(5*s**2 + 54*s*t - 3*t**2)
                + cn.mp**2*(s**3 + 11*s**2*t + 23*s*t**2 + t**3)
                - s*t*(s+t)**2)
        num4 = -2*(4*cn.mp**6 + cn.mp**4*(s-7*t) - 2*cn.mp**2*(s**2-t**2)
                   + s*t*(s+t))
        num6 = 2*(cn.mp**4 - cn.mp**2*(s+t) + s*t)
        return (num0 + num2*mv**2 + num4*mv**4 + num6*mv**6)*denom

    if frame == "com":
        ep = (s + cn.mp**2)/2/np.sqrt(s)
        t = (mv**2 + cn.mp**2 - 2*ev*ep
             - np.sqrt(ev**2-mv**2)*np.sqrt(ep**2-cn.mp**2)*ct)
        return faca * trc(t)
    elif frame == "lab":
        t0 = cn.mp**2 + mv**2 - 2*cn.mp*ev[0]
        t1 = cn.mp**2 + mv**2 - 2*cn.mp*ev[1]
        return (faca * trc(t0), faca * trc(t1))
