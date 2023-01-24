import numpy as np

Eph_lo = 6.
Eph_hi = 24.

rmin = 0.15
rmax = 0.55

zmin = -0.5
zmax = 1.0
Delz_lo = -0.15
Delz_hi = 0.15

eff_el = 0.93
eff_mu = 0.86
eff_hd = 1.

pTlo = 0.9
etlo = -1
ethi = 2.5


def eff(r, z):
    if r > rmin and r < rmax and z > zmin and z < zmax:
        return (rmax - r)/(rmax - rmin)


eta_narrow = np.linspace(etlo, ethi, 500)
