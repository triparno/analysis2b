'''Script to generate differential cross-section graphs'''
import numpy as np
import matplotlib.pyplot as plt
from mat import *
import amp as amp
import matplotlib.ticker as ticker
import cons as cn


def gam_fr(x): return (x+cn.mp)/np.sqrt(cn.mp**2+2*cn.mp*x)


etmin = -np.log(np.tan(20.*np.pi/180/2))
etmax = -np.log(np.tan(130.*np.pi/180/2))
theta_range = np.linspace(0., 180., 10)*np.pi/180
ct_range = np.linspace(-1, 1, 1000)
eta_range = -np.log(np.tan(theta_range/2))

mv_lis = np.logspace(1, 3, 4)/1000.
Egamma_lis = np.arange(6, 13, 6)



ls_lis = [":", "-.", "-", "--"]
lc_lis = ["blue", "red"]
lab_lis = [r"6 GeV", r"12 GeV"]
mas_lis = [r"$10^0$ MeV", r"$10^1$ MeV", r"$10^2$ MeV", r"$10^3$ MeV"]


fig_dir = "fig/"
pp = PdfPages(fig_dir+"DsgDet.pdf")
fig = plt.figure(figsize=(3, 2))
ax = fig.add_subplot()
ax.set_xlabel(r'$\eta$')
ax.set_ylabel(r'$d\sigma/d\eta$')
# ax.set_ylim(0, 100)
ax.set_xlim(-2.5, 2.5)
ii = 0
#   difsec_range_com = amp.AmpSq(np.cos(theta_range), 0.1, 6, 'com')
#   difsec_range_lab = amp.AmpSq(np.cos(theta_range), 0.1, 6, 'lab')
#   plt.plot(theta_range, difsec_range_com, lw=1)
#   plt.plot(theta_range, difsec_range_lab, lw=1)
for ep in Egamma_lis:
    eta_range_translation = eta_range + np.log(gam_fr(ep))
    jj = 0
    for mv in mv_lis:
        difsec_range = amp.dSidet(eta_range, mv, ep, 'com')
        plt.plot(eta_range, difsec_range,
                 ls=ls_lis[jj], c=lc_lis[ii], lw=1)
        difsec_range = amp.dSidet(eta_range, mv, ep, 'lab')
        plt.plot(eta_range, difsec_range,
                 ls=ls_lis[jj], c='k', lw=1)
        jj += 1
    ii += 1

#   for ii in range(0, len(Egamma_lis)):
#       plt.plot(-6, -6.1, lw=1, c=lc_lis[ii], ls='-',
#                label=lab_lis[ii])
#   for ii in range(0, len(mv_lis)):
#       plt.plot(-6, -6.1, lw=1, c='k', ls=ls_lis[ii],
#                label=mas_lis[ii])

#   plt.axvline(etmin, c='0.5', lw=0.5)
#   plt.axvline(etmax, c='0.5', lw=0.5)

plt.yscale('log')
ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=7))

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

# ax.ticklabel_format(axis='y', style='scientific', scilimits=[0, 2])
plt.legend(loc='lower left', fontsize=6)
fig.tight_layout(pad=0.1)
fig.savefig(pp, format='pdf')
pp.close()


mv_lis_dense = np.logspace(-3, 1, 1000)
lenM = len(mv_lis_dense)

totsec_lab = np.zeros(lenM)
for ii in range(0, lenM):
    totsec_lab[ii] = np.trapz(amp.dSidct(ct_range, mv_lis_dense[ii], 12, 'lab'), ct_range)

totsec_com = np.zeros(lenM)
for ii in range(0, lenM):
    totsec_com[ii] = np.trapz(amp.dSidct(ct_range, mv_lis_dense[ii], 12, 'com'), ct_range)

pp = PdfPages(fig_dir+"Sig_mas.pdf")
fig = plt.figure(figsize=(2.7, 2))
ax = fig.add_subplot()
ax.set_xlabel(r'$m_v$ [GeV]')
ax.set_ylabel(r'$\sigma$ [pb]')

# ax.set_ylim(0, 100)
# ax.set_xlim(-2.5, 2.5)

plt.plot(mv_lis_dense, totsec_lab,
         lw=0.5,
         label='lab')
plt.plot(mv_lis_dense, totsec_com,
         lw=0.5,
         label='com')

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')

plt.yscale('log')
plt.xscale('log')

ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))
ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))

# ax.ticklabel_format(axis='y', style='scientific', scilimits=[0, 2])

plt.legend(loc='lower left', fontsize=6)
fig.tight_layout(pad=0.1)
fig.savefig(pp, format='pdf')
pp.close()
