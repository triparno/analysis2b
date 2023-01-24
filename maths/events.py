import numpy as np
import csec as csc
import functions as func
import parm as pr
import cons as cn
from scipy.stats import reciprocal
import inputs

nlo = 50
nin = 900
nsmall = 500

rmin = pr.rmin
rmax = pr.rmax
zmin = pr.zmin
zmax = pr.zmax


def modp(mv, Eph):
    '''Return the magnitude of the three momentum for an exotic vector mass

    Parameters
    ----------
    mv : float
        Mass of the exotic vector
    Eph : float
        Incoming photon Energy
    Returns
    -------
    float
        three momentum magnitude
    '''
    s = cn.mp**2 + 2*Eph*cn.mp
    return np.sqrt(func.lamK(s, mv**2, cn.mp**2)/4/s)


def declen(mv, gg, totgam, Eph):
    '''Return decay length given mass and coupling
    Parameters
    ----------
    mv : float
        Mass of the exotic vector
    gg : float
        Coupling strength of the exotic vector
    totgam : float
        Total width of the exotic vector
    Eph : float
        Energy of the incoming photon
    Returns
    -------
    float
        decay length of the exotic vector
    '''
    ll = cn.cc * cn.GeVIn2s / totgam / gg**2
    return ll * modp(mv, Eph) / mv


def efi_nom(rad, zed):
    '''Return detector efficiency for a particular r and z
    Parameters
    ----------
    rad : float
        radial distance
    zed : float
        distance parallel to beam
    Returns
    -------
    int
        efi_nom, binary 0 or 1
    '''
    return (np.heaviside(zmax-zed, 0)
            * np.heaviside(zed-zmin, 0)
            * np.heaviside(rad-rmin, 0)
            * np.heaviside(rmax-rad, 0)
            )


def eff_ang(mv, gg, totgam, eta, Eph):
    '''Convolute the decay length distribution with the detector
    efficiency distribution
    Parameters
    ----------
    mv : float
        Mass of the exotic vector
    gg : float
        Coupling strength of the exotic vector
    totgam : float
        Total width of the exotic vector
    eta : float
        pseudorapidities
    Eph : float
        Energy of the incoming photon
    Returns
    -------
    float
        eff_ang, decay length distribution convoluted with detector
        efficiency
    '''
    lhalf = declen(mv, gg, totgam, Eph)
    th = 2 * np.arctan(np.exp(-eta))
    lran = np.linspace(rmin, np.hypot(rmin, zmax), 100)
    ldist = np.exp(-lran/lhalf)/lhalf
    norm = np.sum(ldist)
    rr = lran*np.sin(th)
    zz = lran*np.cos(th)
    zz = np.random.uniform(zz - pr.Delz_lo, zz + pr.Delz_hi)
    eff = np.sum(efi_nom(rr, zz) * ldist)/norm/(pr.Delz_hi-pr.Delz_lo)
    return eff


eff_ang_vec = np.vectorize(eff_ang)     # Vectorize the above function


def eta_int(mv, gg, totgam, Ep, frame, eta=pr.eta_narrow):
    '''Perform the eta integral convoluted with efficiency.

    Parameters
    ----------
    mv : float
        Mass of the Dark Photon
    gg : float
        Coupling of the dark photon
    Ep : float
        Incoming photon energy
    eta : numpy array, optional
        Range of pseudorapidity covered (default as given in parameter file)

    Returns
    -------
    float
        eta_int
    '''
    eff = eff_ang_vec(mv, gg, totgam, eta, Ep)
    Dsig = eff * gg**2 * csc.dSidet(eta, mv, Ep, frame)
    Dsig *= np.heaviside(modp(mv, Ep)*np.sin(func.thFet(eta))-pr.pTlo, 0)
    return np.trapz(Dsig, eta)


def select(mv, gg, totgam, frame, BR=1, eff_g=1, lumi=1, Nev=3,
           eta=pr.eta_narrow):
    '''For a particular mass and coupling pair check if number of
    accepted events is greater than required sensitivity. The photon
    energy is selected from a reciprocal distribution.

    Parameters
    ----------
    mv : float
        Mass of the Dark Photon
    gg : float
        Coupling of the dark photon
    BR : float, optional
        Branching Ratio (default is 1)
    eff_g : float, optional
        Overall Efficiency (default is 1)
    lumi : float, optional
        luminosity (default is 1 pbIn)
    Nev : float, optional
        Number of events to cross (default is 1)
    eta : numpy array, optional
        Range of pseudorapidity covered (default as given in parameter file)

    Returns
    -------
    int
        binary, 1 or zero
    '''
    arrE = np.linspace(pr.Eph_lo, pr.Eph_hi, inputs.Eph_st)
    pdfE = reciprocal.pdf(arrE, pr.Eph_lo, pr.Eph_hi)
    pdfE /= np.sum(pdfE)
    Eph = np.random.choice(arrE, p=pdfE)
    tot = eta_int(mv, gg, totgam, Eph, frame)
    tot *= lumi * BR * eff_g
    return np.heaviside(tot-Nev, 0)
