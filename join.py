import numpy as np
import inputs

model = inputs.model
ch = inputs.ch
lum =inputs.lum
Nevents = inputs.Nevents

dathi1 = np.genfromtxt("data/"+model+"_"+ch+f"_{lum}_{Nevents}_hi_tp.dat").T
datlo1 = np.genfromtxt("data/"+model+"_"+ch+f"_{lum}_{Nevents}_lo_tp.dat").T
# dathi2 = np.genfromtxt(model+"_"+ch+f"_{lum}_{Nevents}_hi_cl.dat").T
# datlo2 = np.genfromtxt(model+"_"+ch+f"_{lum}_{Nevents}_lo_cl.dat").T


with open(model+"_"+ch+f"_{lum}_{Nevents}.csv", "w") as out:
    for ii in range(len(datlo1[0])):
        out.write(f"\n{datlo1[0][ii]} {datlo1[1][ii]} {dathi1[1][ii]}")
