from modules import dat_red
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from had_dec import HDout
from matplotlib.colors import LogNorm
from mat import *

had_path = "/home/triparno/Dropbox/b2sz/code/data/BR_full/had/"     # Location of hadron data
tau_path = "/home/triparno/Dropbox/b2sz/code/data/BR_full/tau/"     # Location of tau data
elc_path = "/home/triparno/Dropbox/b2sz/code/data/BR_full/elc/"     # Location of electron data
muo_path = "/home/triparno/Dropbox/b2sz/code/data/BR_full/muo/"     # Location of muon data
DC_path = "/home/triparno/Dropbox/b2sz/code/data/BR_full/DC/"       # Location of DC BF outputs

DC_ll_BF = np.genfromtxt(DC_path+"bfrac_dark_photon_ell.txt").T
DC_hh_BF = np.genfromtxt(DC_path+"bfrac_dark_photon_hadrons.txt").T
DC_tt_BF = np.genfromtxt(DC_path+"bfrac_dark_photon_tau_tau.txt").T

clp_x = [0., 25.]
clp_y = [0., 1.]

ah, bh = dat_red.datTrcReg(had_path, clp_x, clp_y)
ae, be = dat_red.datTrcReg(elc_path, clp_x, clp_y)
am, bm = dat_red.datTrcReg(muo_path, clp_x, clp_y)
at, bt = dat_red.datTrcReg(tau_path, clp_x, clp_y)

GeVIn2s = 6.58*10**(-25)
cc = 3*10**8
aEM = 1/137.
me = 0.5E-3
mm = 1.0E-3
mB = 5.3
mK = 4.9E-3
eEM = np.sqrt(4*np.pi*aEM)


def lamK(a,b,c):
    return np.sqrt(a**4+b**4+c**4-2*(a**2*b**2+b**2*c**2+c**2*a**2))


def D2fer(mX,ep,ml,nc):
    '''
    Decay Rate into fermions
    '''
    cc= np.heaviside(mX-2*ml,0)*(1-4*ml**2/mX**2)
    return aEM*nc*ep**2*mX*np.sqrt(cc)*(1+2*ml**2/mX**2)/3


mtu=1.7
mbt= 4.2
mch= 1.3
mst= 0.95

Dtu= lambda mX,ep: D2fer(mX,ep,mtu,1) 
Dmu= lambda mX,ep: D2fer(mX,ep,0.1,1)
De= lambda mX,ep: D2fer(mX,ep,0.,1)

Dbt= lambda mX,ep: D2fer(mX,ep,4.2,3)/9
Dch= lambda mX,ep: D2fer(mX,ep,mch,3)*4/9
Dst= lambda mX,ep: D2fer(mX,ep,mst,3)/9
Dup= lambda mX,ep: D2fer(mX,ep,0.,3)*4/9
Ddn= lambda mX,ep: D2fer(mX,ep,0.,3)/9



def D2lepT(mX,ep):
    return Dmu(mX,ep)+De(mX,ep)

def D2hadT(mX,ep):
    return Dup(mX,ep)+Ddn(mX,ep)+Dbt(mX,ep)+Dch(mX,ep)+Dst(mX,ep)

def D2allT(mX,ep):
    return Dtu(mX,ep)+Dmu(mX,ep)+De(mX,ep)+Dbt(mX,ep)+Dch(mX,ep)+Dst(mX,ep)+Dup(mX,ep)+Ddn(mX,ep)


#def D2had(mX,g):
#    '''
#    Decay Rate into hadrons
#    '''
#    return (10**g)**2*mX*HDout(mX)/12/np.pi

#def D2all(mX,g):
#    '''
#    Total Decay Rate 
#    '''
#    return D2lep(mX,me,g)+D2lep(mX,mm,g)+D2had(mX,g)

def BRtu_pQCD(mX,g):
    return Dtu(mX,g)/D2allT(mX,g)

def BRlep_pQCD(mX,g):
    return D2lepT(mX,g)/D2allT(mX,g)

def BRhad_pQCD(mX,g):
    return D2hadT(mX,g)/D2allT(mX,g)

def BRDL(mX,g):
    '''
    Branching to leptons 
    '''
    return (D2lep(mX,me,g)+D2lep(mX,mm,g))/D2all(mX,g)

def BRDH(mX,g):
    '''
    Branching to leptons 
    '''
    return D2had(mX,g)/D2all(mX,g)

def DL(mX,g):
    sqlamK= np.sqrt(mB**4+mK**4+mX**4-2*mB**2*mK**2-2*mB**2*mX**2-2*mK**2*mX**2)
    mp= sqlamK/2/mB; btgm= mp/mX
    return btgm*GeVIn2s*cc/D2all(mX,g)

def out(mX,gx):
    NN= (55/2.)*1E9
    BKga= 3.3E-4
    #*lamK(mB,mK,mX)/lamK(mB,mK,0)
    #return (NN)*BRDL(mX,gx)*(BKga*(10**gx)**2/(4*np.pi*aEM))
    return (BKga*(10**gx)/(4*np.pi*aEM))


mm=(min(ae),min(am))
xx=(max(ae),max(am))

ll= (max(mm),min(xx))

ranl= np.linspace(ll[0],ll[1],100)

inte= interpolate.interp1d(ae,be)
intm= interpolate.interp1d(am,bm)

mout= intm(ranl)
eout= inte(ranl)

lout= mout + eout 
#print(mout)

pp=PdfPages('BR_pqcd.pdf')
fig=plt.figure(figsize=(4,3))
ax= plt.axes()

mvar= np.linspace(0.3,20,100)


plt.xlim(0.0,12)

plt.plot(mvar,BRhad_pQCD(mvar,1),ls='--',c=afblue)
plt.plot(ah,bh, c=afblue,label='hadron')
plt.plot(DC_hh_BF[0], DC_hh_BF[1],c=afblue,ls=':')


plt.plot(ranl,lout,label='lepton', c=sig1)
plt.plot(mvar,BRlep_pQCD(mvar,1), c=sig1,ls='--')
plt.plot(DC_ll_BF[0], DC_ll_BF[1],c=sig1,ls=':')

plt.plot(mvar,BRtu_pQCD(mvar,1),c=red,ls='--')
plt.plot(at,bt,label='tau',c=red)
plt.plot(DC_tt_BF[0], DC_tt_BF[1],c=red,ls=':')

plt.legend()
ax.set_xlabel(r'$M_{\gamma^\prime}$ [GeV]')
ax.set_ylabel(r'BR')
#ax.set_ylabel(r'$\log_{10}(\epsilon e)$')
#fig.colorbar(ll,label=r'Decay Length [m]')

fig.tight_layout()
fig.savefig(pp,format='pdf')
pp.close()
