# ------------------------------------------------------------------- #
# ----------- Script to read the traced data corresponding  --------- #
# ----------- to the branching ratios of the dark photon    --------- #
# ----------- as given in arxiv:1801.04847 [hep-ph].        --------- #
# ------------------------------------------------------------------- #
from modules import mat
import pandas as pd
import os
from matplotlib import pyplot as plt

data_path="/home/triparno/Dropbox/b2sz/code/data/"

# We will read three different traces corresponding to the three 
# channels and average over them (for each channel)
ee,mm,hh= [],[],[]
for ii in (1,2,3):
    ee.append(pd.read_csv(data_path+"SDP/ee%i.csv"%(ii))) # el channel
    mm.append(pd.read_csv(data_path+"SDP/mm%i.csv"%(ii))) # mu channel
    hh.append(pd.read_csv(data_path+"SDP/hh%i.csv"%(ii))) # hd channel

# Trimming the data to handle points outside the axes
cl, cr, cb, ct= 0,2,0,1
for ii in (0,1,2):
    ee[ii].loc[ee[ii].Curve1>ct,"Curve1"]=ct
    ee[ii].loc[ee[ii].Curve1<cb,"Curve1"]=cb
    ee[ii].loc[ee[ii].x<cl,"x"]=cl
    ee[ii].loc[ee[ii].x>cr,"x"]=cr
    mm[ii].loc[mm[ii].Curve1>ct,"Curve1"]=ct
    mm[ii].loc[mm[ii].Curve1<cb,"Curve1"]=cb
    mm[ii].loc[mm[ii].x<cl,"x"]=cl
    mm[ii].loc[mm[ii].x>cr,"x"]=cr
    hh[ii].loc[hh[ii].Curve1>ct,"Curve1"]=ct
    hh[ii].loc[hh[ii].Curve1<cb,"Curve1"]=cb
    hh[ii].loc[hh[ii].x<cl,"x"]=cl
    hh[ii].loc[hh[ii].x>cr,"x"]=cr


# Arranging data to assure proper x-ordering
ees,mms,hhs= [],[],[]
for ii in (0,1,2):
    ees.append(ee[ii].sort_values(by= 'x', ascending='True'))
    mms.append(mm[ii].sort_values(by= 'x', ascending='True'))
    hhs.append(hh[ii].sort_values(by= 'x', ascending='True'))


# Getting the proper endpoints for consistent interpolation
eehi, eelo, mmhi, mmlo, hhhi, hhlo= [],[],[],[],[],[]
for ii in (0,1,2):
    eelo.append(min(ees[ii].x)); eehi.append(max(ees[ii].x))
    mmlo.append(min(mms[ii].x)); mmhi.append(max(mms[ii].x))
    hhlo.append(min(hhs[ii].x)); hhhi.append(max(hhs[ii].x))

elim= [max(eelo),min(eehi)]; mlim= [max(mmlo),min(mmhi)]; hlim= [max(hhlo),min(hhhi)]

ii=0
while ii<2:
    del ee[ii]; del mm[ii]; del hh[ii]
    ii+= 1

# Interpolation between endpoints obtained in last step 
eei,mmi,hhi= [],[],[]
for ii in (0,1,2):
    eei.append(interpolate.interp1d(ees[ii].x,ees[ii].Curve1))
    mmi.append(interpolate.interp1d(mms[ii].x,mms[ii].Curve1))
    hhi.append(interpolate.interp1d(hhs[ii].x,hhs[ii].Curve1))

eelen= np.linspace(elim[0],elim[-1],100)
mmlen= np.linspace(mlim[0],mlim[-1],100)
hhlen= np.linspace(hlim[0],hlim[-1],100)

# Getting a final dataset after averaging over the three traces 
# corresponding to the three channels 
eeht, mmht, hhht= [], [], []
for i in range(0,len(eelen)):
    eeht.append((eei[0](eelen[i])+eei[1](eelen[i])+eei[2](eelen[i]))/3.)
    mmht.append((mmi[0](mmlen[i])+mmi[1](mmlen[i])+mmi[2](mmlen[i]))/3.)
    hhht.append((hhi[0](hhlen[i])+hhi[1](hhlen[i])+hhi[2](hhlen[i]))/3.)

