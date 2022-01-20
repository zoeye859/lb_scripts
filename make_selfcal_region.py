from __future__ import division
from astropy import units as u
from astropy.coordinates import SkyCoord
import os
import numpy as np
import re
from regions import DS9Parser, write_ds9

def filename_to_skycoord(filename):
  [left, right] = re.findall(r'(\d+.\d+)', filename)
  hh = left[0:2]
  mm = left[2:4]
  ss = left[4:]
  dd = right[0:2]  
  m2 = right[2:4]
  s2 = right[4:]
  region_centre = '{}h{}m{}s +{}d{}m{}s'.format(hh, mm, ss, dd, m2, s2)
  phasecenter = SkyCoord(region_centre, frame='fk5')
  
  return phasecenter.ra, phasecenter.dec

file_path = '/net/rijn3/data2/sweijen/en1/ddecals/'
files = os.listdir(file_path)
region_list = []

for file in files:
  if '.ms' in file:
    [RA, DEC] = filename_to_skycoord(file)
    region_list.append((RA.deg,DEC.deg))
    region_strs = list(map(lambda pos: 'fk5\npoint({:f},{:f})# color=red width=1 text=""'.format(*pos), region_list))

   
parser = DS9Parser('\n'.join(list(region_strs)))
regions = parser.shapes.to_regions() 
write_ds9(regions, 'selfcal_regions.reg')
    