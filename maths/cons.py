import numpy as np

GeVIn2s= 6.58*10**(-25)
GeVInSq2pb= 0.39*1E9
cc= 3*10**8     # m/s
aEM= 1/137.     # Fine Structure constant at 
me= 0.5E-3      # GeV Electron Mass
mm= 1E-1      # GeV Muon Mass
mt= 1.8         # GeV Muon Mass
mB= 5.3         # GeV B0   Mass
mK= 4.9E-3      # GeV K    Mass
mp = 0.98 
eEM= np.sqrt(4*np.pi*aEM)
NN= (55/2.)*1E9     # Number of bbbar collisions at Belle II
BKga= 3.3E-4        # Branching B to X gamma

def lamK(a,b,c):
    return a**2+b**2+c**2-2*(a*b+b*c+c*a)
