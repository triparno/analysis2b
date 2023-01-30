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

data_fp_u = np.genfromtxt(inputs.cnt_file+"_sp_u.dat").T
data_fp_d = np.genfromtxt(inputs.cnt_file+"_sp_d.dat").T

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

steps_u = range(len(data_fp_u[0])-1)
steps_d = range(len(data_fp_d[0])-1)
mlis_u = data_fp_u[0]
mlis_d = data_fp_d[0]

gstep = (np.log10(inputs.grange[1])-np.log10(inputs.grange[0]))/10


def ff_u(ii):
    lismu = []
    lisgu = []
    gval = np.log10(data_fp_u[1][ii])
    gvalmin, gvalmax = gval, gval + gstep
    mm = data_fp_u[0][ii]
    br = br_func(ch, mm)
    for gg in np.logspace(gvalmin, gvalmax, inputs.gscan_step):
        if evn.select(mm, inputs.coup_scal*gg, IF[0](mm), pr.eta_narrow,
                      frame, BR=br, eff_g=global_ef, lumi=inputs.lum,
                      Nev=inputs.Nevents) == 0:
            lismu.append(mm)
            lisgu.append(gg)
            break
    return (lismu, lisgu)


def ff_d(ii):
    lismd = []
    lisgd = []
    gval = np.log10(data_fp_d[1][ii])
    gvalmin, gvalmax = gval - gstep, gval
    mm = data_fp_d[0][ii]
    br = br_func(ch, mm)
    for gg in np.logspace(gvalmin, gvalmax, inputs.gscan_step):
        if evn.select(mm, inputs.coup_scal*gg, IF[0](mm), pr.eta_narrow,
                      frame, BR=br, eff_g=global_ef, lumi=inputs.lum,
                      Nev=inputs.Nevents) == 1:
            lismd.append(mm)
            lisgd.append(gg)
            break
    return (lismd, lisgd)


if __name__ == '__main__':
    with multiprocessing.Pool(12) as p:
        lis_u = p.map(ff_u, steps_u)

if __name__ == '__main__':
    with multiprocessing.Pool(12) as p:
        lis_d = p.map(ff_d, steps_d)

with open(inputs.cnt_file+"_tp_u.dat", "w") as out:
    for ii in range(len(lis_u)):
        for jj in range(len(lis_u[ii][0])):
            out.write(f"\n{lis_u[ii][0][jj]} {lis_u[ii][1][jj]}")

with open(inputs.cnt_file+"_tp_d.dat", "w") as out:
    for ii in range(len(lis_d)):
        for jj in range(len(lis_d[ii][0])):
            out.write(f"\n{lis_d[ii][0][jj]} {lis_d[ii][1][jj]}")
et = time.time()
print(f"Total time is: {et-st}")
