import numpy as np
import parm as pr
import events as evn
import decay as dec
import inputs

ch = inputs.ch

IF = dec.Gam_get(inputs.model, datW=inputs.wd_dir)

data_fp = np.genfromtxt(inputs.cnt_file+"_sp.dat").T

frame = inputs.frame


def br_func(ch, mv):
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
gstep = np.log10(inputs.grange[1])-np.log10(inputs.grange[0])

mLislo, mLishi, gLislo, gLishi = [], [], [], []

for ii in range(len(data_fp[0])):
    gval = np.log10(data_fp[2][ii])
    mm = data_fp[0][ii]
    br = br_func(ch, mm)
    gvalmin, gvalmax = gval, gval + gstep
    for gg in np.logspace(gvalmax, gvalmin, inputs.scan_step):
        if evn.select(mm, inputs.coup_scal*gg, IF[0](mm), pr.eta_narrow, frame, BR=br, eff_g=global_ef, lumi=inputs.lum,
                      Nev=inputs.Nevents) == 1:
            mLishi.append(mm)
            gLishi.append(gg)
            break

for ii in range(len(data_fp[0])):
    gval = np.log10(data_fp[1][ii])
    mm = data_fp[0][ii]
    br = br_func(ch, mm)
    gvalmin, gvalmax = gval - gstep, gval
    for gg in np.logspace(gvalmin, gvalmax, inputs.scan_step):
        if evn.select(mm, inputs.coup_scal*gg, IF[0](mm), pr.eta_narrow, frame, BR=br, eff_g=global_ef, lumi=inputs.lum,
                      Nev=inputs.Nevents) == 1:
            mLislo.append(mm)
            gLislo.append(gg)
            break

with open(inputs.cnt_file+"_tp.dat", "w") as out:
    for ii in range(len(mLislo)):
        out.write(f"\n{mLislo[ii]} {gLislo[ii]} {gLishi[ii]}")
