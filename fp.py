import time
import numpy as np
import decay as dec
import parm as pr
import events as evn
import inputs


st = time.time()
ch = inputs.ch
frame = inputs.frame

IF = dec.Gam_get(inputs.model, datW=inputs.wd_dir)


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
elif ch == 'mu':
    global_ef = pr.eff_mu
elif ch == 'hd':
    global_ef = pr.eff_hd
elif ch == 'lp':
    global_ef = pr.eff_mu

mLislo, mLishi, gLislo, gLishi = [], [], [], []

for mm in inputs.mrange:
    br = br_func(ch, mm)
    for gg in np.flipud(inputs.grange):
        cp = inputs.coup_scal*gg
        csec = evn.select(mm, cp, IF[0](mm), pr.eta_narrow, frame, BR=br, eff_g=global_ef,
                          lumi=inputs.lum, Nev=inputs.Nevents)
        if csec == 1:
            mLishi.append(mm)
            gLishi.append(gg)
            break

for mm in inputs.mrange:
    br = br_func(ch, mm)
    for gg in inputs.grange:
        cp = inputs.coup_scal*gg
        csec = evn.select(mm, cp, IF[0](mm), pr.eta_narrow, frame, BR=br, eff_g=global_ef,
                          lumi=inputs.lum, Nev=inputs.Nevents)
        if csec == 1:
            mLislo.append(mm)
            gLislo.append(gg)
            break

with open(inputs.cnt_file+"_fp.dat", "w") as out:
    for ii in range(len(mLislo)):
        out.write(f"\n{mLislo[ii]} {gLislo[ii]} {gLishi[ii]}")
et = time.time()
print(f"Total time is: {et-st}")
