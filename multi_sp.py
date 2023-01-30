import numpy as np
import parm as pr
import events as evn
import decay as dec
import inputs
import multiprocessing
import time
st = time.time()

ch = inputs.ch

IF = dec.Gam_get(inputs.model, datW=inputs.wd_dir)

data_fp = np.genfromtxt(inputs.cnt_file+"_fp.dat").T

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

steps = range(len(data_fp[0])-1)
mlis = data_fp[0]


def ff(ii):
    lisgu = []
    lismu = []
    lisgd = []
    lismd = []
    gval = np.log10(data_fp[1][ii])
    gvalmin, gvalmax = gval, gval + gstep
    mm = data_fp[0][ii]
    br = br_func(ch, mm)
    for gg in np.logspace(gvalmin, gvalmax, inputs.gscan_step):
        if evn.select(mm, inputs.coup_scal*gg, IF[0](mm), pr.eta_narrow,
                      frame, BR=br, eff_g=global_ef, lumi=inputs.lum,
                      Nev=inputs.Nevents) == 0:
            lisgu.append(gg)
            lismu.append(mm)
            break
    gval = np.log10(data_fp[2][ii])
    gvalmin, gvalmax = gval - gstep, gval
    for gg in np.logspace(gvalmin, gvalmax, inputs.gscan_step):
        if evn.select(mm, inputs.coup_scal*gg, IF[0](mm), pr.eta_narrow,
                      frame, BR=br, eff_g=global_ef, lumi=inputs.lum,
                      Nev=inputs.Nevents) == 1:
            lisgd.append(gg)
            lismd.append(mm)
            break
    return (lismu, lisgu, lismd, lisgd)


if __name__ == '__main__':
    with multiprocessing.Pool(12) as p:
        lis = p.map(ff, steps)

with open(inputs.cnt_file+"_sp_u.dat", "w") as out:
    for ii in range(len(lis)):
        for jj in range(len(lis[ii][0])):
            out.write(f"\n{lis[ii][0][jj]} {lis[ii][1][jj]}")

with open(inputs.cnt_file+"_sp_d.dat", "w") as out:
    for ii in range(len(lis)):
        for jj in range(len(lis[ii][2])):
            out.write(f"\n{lis[ii][2][jj]} {lis[ii][3][jj]}")
et = time.time()
print(f"Total time is: {et-st}")
