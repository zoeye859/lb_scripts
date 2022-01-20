from __future__ import division
from astropy import units as u
from astropy.coordinates import SkyCoord
import pyregion
from astropy.io import fits
import os

import numpy as np

phasecenter = SkyCoord('16h11m00s +54d57m00s', frame='fk5')
Phase_centre_RA = (360 - 117.25) * u.degree
spacing = 0.5 * u.degree
# Have some overlap
facet_size = 0.55 * u.degree
facets = np.zeros((5,5,2)) * u.degree
facetlist = []
facetlist_uncorr = []
k = 1
for i in range(5):
    for j in range(5):
        RA = phasecenter.ra + (spacing * (j-2.) / np.cos((phasecenter.dec + spacing*(i-2)).rad))
        RA2 = phasecenter.ra + (spacing * (j-2.))
        DEC = phasecenter.dec + (spacing * (i-2.))
        angle = phasecenter.ra - RA
        facets[i, j, 0] = RA
        facets[i, j, 1] = DEC
        facetlist.append((RA.deg,DEC.deg,facet_size.value,facet_size.value, angle.deg))
        PARSET = 'msout.storagemanager=dysco\nmsout.storagemanager.databitrate=4\nmsout.storagemanager.weightbitrate=8\nsteps=[shift,avg]\nshift.type = phaseshift\nshift.phasecenter = [{:f}deg, {:f}deg]\navg.type = average\navg.timeresolution = 4\navg.freqresolution = 48.82kHz'.format(RA.deg, DEC.deg)
        with open('shift_to_facet_{:d}.parset'.format(k), 'w') as f:
            f.write(PARSET)
        k += 1
#region_strs = list(map(lambda pos: 'fk5\nbox({:f},{:f},{:f},{:f},0) # color=green width=4 text=""'.format(*pos, facet_size.value, facet_size.value), facetlist)) # this only works for python3
region_strs = list(map(lambda pos: 'fk5\nbox({:f},{:f},{:f},{:f},{:f}) # color=green width=4 text=""'.format(*pos), facetlist))

for i, f in enumerate(region_strs):
    with open('facet_{:02d}.reg'.format(i+1), 'w') as ff:
        ff.write(f)

# Make a single region file for easier overlaying in ds9.
from regions import DS9Parser, write_ds9
parser = DS9Parser('\n'.join(list(region_strs)))
regions = parser.shapes.to_regions()

if os.path.isfile('facets.reg'):
      os.system('rm -rf facets.reg') 
      
write_ds9(regions, 'facets.reg')

