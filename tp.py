import numpy as np
import amp as amp
import parm as pr
import events as evn
import inputs

model = inputs.model
IF = amp.Gam_get(model)
ch = inputs.ch
lum = inputs.lum
Nevents = inputs.Nevents

IF = amp.Gam_get(model)

dathi = np.genfromtxt("data/"+model+"_"+ch+f"_{lum}_{Nevents}_hi_sp.dat").T
datlo = np.genfromtxt("data/"+model+"_"+ch+f"_{lum}_{Nevents}_lo_sp.dat").T


def inputs(ch, mv):
    if ch == 'el':
        BR = IF[1](mv)/IF[0](mv)
    if ch == 'mu':
        BR = IF[2](mv)/IF[0](mv)
    if ch == 'hd':
        BR = IF[-1](mv)/IF[0](mv)
    if ch == 'lp':
        BR = (IF[1](mv)+IF[2](mv))/IF[0](mv)
    return BR


if ch == 'el':
    global_ef = pr.eff_el
if ch == 'mu':
    global_ef = pr.eff_mu
if ch == 'hd':
    global_ef = pr.eff_hd
if ch == 'lp':
    global_ef = pr.eff_mu

l10 = np.log(10.)


with open("data/"+model+"_"+ch+f"_{lum}_{Nevents}_hi_tp.dat", "w") as out:
    for ii in range(len(dathi[0])):
        gval = np.log(dathi[1][ii])
        mm = dathi[0][ii]
        br = inputs(ch, mm)
        for gg in np.logspace(gval/l10+0.25, gval/l10, 50):
            if evn.select(mm, gg, BR=br, eff_g=global_ef, lumi=lum, Nev=Nevents) == 1:
                out.write(f"\n{mm} {gg}")
                break
with open("data/"+model+"_"+ch+f"_{lum}_{Nevents}_lo_tp.dat", "w") as out:
    for ii in range(len(datlo[0])):
        gval = np.log(datlo[1][ii])
        mm = datlo[0][ii]
        br = inputs(ch, mm)
        for gg in np.logspace(gval/l10-0.25, gval/l10, 50):
            if evn.select(mm, gg, BR=br, eff_g=global_ef, lumi=lum, Nev=Nevents) == 1:
                out.write(f"\n{mm} {gg}")
                break
