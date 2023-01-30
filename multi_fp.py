import time
import numpy as np
import decay as dec
import parm as pr
import events as evn
import inputs
import multiprocessing


st = time.time()
ch = inputs.ch
frame = inputs.frame
mlis = inputs.mrange

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


def ff(mm):
    br = br_func(ch, mm)
    csec = 0
    gg1 = 10
    gg2 = 10
    for ggu in np.flipud(inputs.grange):
        cp = inputs.coup_scal*ggu
        csec = evn.select(mm, cp, IF[0](mm), pr.eta_narrow, frame,
                          BR=br, eff_g=global_ef, lumi=inputs.lum,
                          Nev=inputs.Nevents)
        if csec == 1:
            gg1 = ggu
            break
    csec = 0
    for ggd in inputs.grange:
        cp = inputs.coup_scal*ggd
        csec = evn.select(mm, cp, IF[0](mm), pr.eta_narrow, frame,
                          BR=br, eff_g=global_ef, lumi=inputs.lum,
                          Nev=inputs.Nevents)
        if csec == 1:
            gg2 = ggd
            break
    return (gg1, gg2)


if __name__ == '__main__':
    with multiprocessing.Pool(12) as p:
        lis = p.map(ff, mlis)


with open(inputs.cnt_file+"_fp.dat", "w") as out:
    counter = 0
    for ii in range(len(lis)):
        if lis[ii][0] != 10:
            if lis[ii][1] != 10:
                counter += 1
                out.write(f"\n{mlis[ii]} {lis[ii][0]} {lis[ii][1]}")
    out.write(f"\n{mlis[counter]} {lis[counter-1][0]} {lis[counter-1][1]}")
et = time.time()
print(f"Total time is: {et-st}")
