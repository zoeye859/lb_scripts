#! /usr/bin/env python2

import sys
#import pydal as dal
import matplotlib.pyplot as pl
import pyrap.tables as pt
#import utils.plot_util as pp
import numpy as np
import matplotlib
#import plotting as pp
import os

def get_uvw_freq(MS, freq_samp):
    chantable = pt.table(MS+'/SPECTRAL_WINDOW/')
    ffreq_col = chantable.getcol('CHAN_FREQ')
    ffreq_col = ffreq_col.flatten()
    print "{n} freqeuncies".format(n=len(ffreq_col))
    freq_col = ffreq_col[0:len(ffreq_col):freq_samp]
    return freq_col


def get_uvw(MS, time_samp):
    maintable = pt.table(MS)
    fuvw_column = maintable.getcol("UVW")
    print "{n} times".format(n=len(fuvw_column))
    uvw_column = fuvw_column[0:len(fuvw_column):time_samp,:]
    vdata=uvw_column[:,1]
    udata=uvw_column[:,0]
    return vdata, udata


directory = '/net/rijn2/data2/Haoyang/Elais/backup/ELAIS-N1/'
MS = '/net/rijn2/data2/Haoyang/Elais/backup/ELAIS-N1/sub6asec_L686962_SB001_uv_12CFFDDA8t_121MHz.ms.sub.shift.avg.apply_infield.avg.ms'
figname = 'uvcoverage'
color_list = pl.cm.BrBG(np.linspace(0, 1, 2))
color_list = pl.cm.cubehelix(np.linspace(0, 1, 4))
cols = color_list[1:]

freq_samp = 40
time_samp = 10 

vdata, udata = get_uvw(MS, time_samp)


freq_col = []
for filename in os.listdir(directory):
    if filename.startswith("sub6asec_L686962_SB001_uv_12CFFDDA8t_") and filename.endswith("MHz.ms.sub.shift.avg.apply_infield.avg.ms"):
        print ("The current .ms file is " + filename + "\n")
        freq_col_temp = get_uvw_freq(directory+filename, time_samp)
        freq_col += freq_col_temp.tolist()
        
Npnts = len(udata)*len(freq_col)
print "selection: {n}".format(n=figname)
print "{n} freqeuncies".format(n=len(freq_col))
print freq_col
print "{n} times".format(n=len(udata))
print "{n} points to be plotted".format(n=Npnts)

    
f=pl.figure()
ax = pl.subplot(111)
#f,ax = pp.paper_single_ax()
pl.minorticks_on()
pl.axis('equal')
# plot the data
ax.set_xlabel("$u$ [k$\lambda$]")
ax.set_ylabel("$v$ [k$\lambda$]")
for freq in freq_col:
    wave = 1000*2.998e8/freq  # freq in Hz
    udatawave = udata/wave
    vdatawave = vdata/wave
    ax.plot(udatawave,vdatawave,'.',c=cols[0], markersize=2, alpha=0.1)
    ax.plot(-udatawave,-vdatawave,'.',c=cols[1], markersize=2, alpha=0.1)
pl.savefig('uv_plot.png', dpi = 300)


