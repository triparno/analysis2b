import numpy as np
import modules.amp as amp
import modules.parm as pr
import modules.cons as cn


model = "dp"   # input("Model name: ")
IF = amp.Gam_get(model)

nlo = 50
nin = 900

eta_lo = np.linspace(-3, -1, nlo)
eta_md = np.linspace(-1, 3, nin)
eta_hi = np.linspace(3, 4, nlo)
eta_ran = np.concatenate((eta_lo, eta_md, eta_hi))


def declen(mv, gg):
    mod_p = (pr.S2-mv**2)/2/np.sqrt(pr.S2)
    ll = cn.cc*cn.GeVIn2s/IF[0](mv)/gg**2
    return ll * mod_p / mv


def efi_nom(rad, zed):
    return (np.heaviside(pr.zmax-zed, 0)
            * np.heaviside(zed-pr.zmin, 0)
            * np.heaviside(rad-pr.rmin, 0)
            * np.heaviside(pr.rmax-rad, 0)
            )  # *(pr.rmax-rad)/(pr.rmax-pr.rmin)


def eff_ang(mv, gg, eta):
    lhalf = declen(mv, gg)
    th = 2*np.arctan(np.exp(-eta))
    lran = np.linspace(pr.rmin, np.hypot(pr.rmin, pr.rmax), 100)
    ldist = np.exp(-lran/lhalf)/lhalf
    norm = np.sum(ldist)
    rr = lran*np.sin(th)
    zz = lran*np.cos(th)
    eff = np.sum(efi_nom(rr, zz)*ldist)
    return eff/norm


eff_ang_vec = np.vectorize(eff_ang)


def eff(mv, gg, eta=eta_ran):
    lhalf = declen(mv, gg)
    th = 2*np.arctan(np.exp(-eta))
    lran = np.linspace(pr.rmin, np.hypot(pr.rmin, pr.rmax), 100)
    ldist = np.exp(-lran/lhalf)/lhalf
    eff = []
    for an in th:
        rr = lran*np.sin(an)
        zz = lran*np.cos(an)
        eff0 = 0
        effSc = 0
        for ii in range(len(rr)):
            eff0 += efi_nom(rr[ii], zz[ii])*ldist[ii]
            effSc += ldist[ii]
        eff.append(eff0/effSc)
    return np.asarray(eff)

def eff_old(mv, gg, eta=eta_ran):
    th = 2*np.arctan(np.exp(-eta))
    mp = (pr.S2-mv**2)/2/np.sqrt(pr.S2)
    ll = cn.cc*cn.GeVIn2s/IF[0](mv)/gg**2
    ll *= mp/mv
    rr = ll*np.sin(th)
    zz = ll*np.cos(th)
    eff0 = (np.heaviside(pr.zmax-zz, 0)
            * np.heaviside(zz-pr.zmin, 0)
            * np.heaviside(rr-pr.rmin, 0)
            * np.heaviside(pr.rmax-rr, 0)
            )*(pr.rmax-rr)/(pr.rmax-pr.rmin)
    return eff0

mm = 0.6
gg = 1e-5
ee = 0.

#   for mm in np.logspace(-2, 0.5, 10):
#       for gg in np.logspace(-6, -3, 10):
#           pp = PdfPages(f'/home/triparno/Projects/b2sz/pt/eff_out/cont_throw_{mm}_{gg}.pdf')
#           fig = plt.figure(figsize=(4, 3))
#           ax = plt.axes()

#           plt.plot(eta_ran, eff_ang_vec(mm, gg, eta_ran), c='blue')
#           plt.plot(eta_ran, eff_old(mm, gg), c='red')
#           plt.axhline(1, c='black')

#           ax.set_xlabel(r'$\eta$')
#           ax.set_ylabel(r'$\mathcal{E}(\eta)$')
#           # ax.set_ylim(1E-7,1E-2)
#           ax.set_xlim(np.min(eta_ran), np.max(eta_ran))
#           plt.text(np.min(eta_ran), 0.9, f'{mm} {gg}')
#           plt.text(np.min(eta_ran), 0.1, f'{declen(mm, gg)}')
#           plt.legend(loc='best')
#           fig.tight_layout()
#           fig.savefig(pp, format='pdf')
#           pp.close()

def select_had(mv, gg, eta=eta_ran):
    BR = IF[-1](mv)/IF[0](mv)
    eff = eff_ang_vec(mv, gg, eta)
    Dsig = eff * amp.dSidetCM(eta, mv, gg) * BR
    tot = np.trapz(Dsig, eta)*50*1.98**2*10**14
    return np.heaviside(tot-3, 0)

def select_lep(mv, gg, eta=eta_ran):
    BR = (IF[1](mv)+IF[2](mv))/IF[0](mv)
    eff = eff_ang_vec(mv, gg, eta)
    Dsig = eff * amp.dSidetCM(eta, mv, gg) * BR
    tot = np.trapz(Dsig, eta)*50*1.98**2*10**14
    return np.heaviside(tot-3, 0)

def select_el(mv, gg, eta=eta_ran):
    BR = IF[1](mv)/IF[0](mv)
    eff = eff_ang_vec(mv, gg, eta)
    Dsig = eff * amp.dSidetCM(eta, mv, gg) * BR * pr.eff_el
    tot = np.trapz(Dsig, eta)*50*1.98**2*10**14
    return np.heaviside(tot-3, 0)

#   start = time.time()
#   with open("pscan_"+model+"_lep_throw_hi_fp.dat", "w") as out:
#       for mm in np.logspace(-2, 0.5, 50):
#           ll = []
#           for gg in np.logspace(-3, -6, 50):
#               if select_lep(mm, gg) == 1:
#                   gmin = gg
#                   out.write(f"\n{mm} {gmin}")
#                   break
#   end = time.time()
#   print(end-start)


with open("pscan_"+model+"_el_throw.dat", "w") as out:
    for mm in np.logspace(-2, 0.4, 50):
        ll = []
        for gg in np.logspace(-2, -7, 50):
            if select_el(mm, gg) == 1:
                ll.append(gg)
        if len(ll) > 0:
            gmin = min(ll)
            gmax = max(ll)
            out.write(f"\n{mm} {gmin} {gmax}")
