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

def num_filename(dir_num, dir_path):
  dir_num = str(dir_num)
  files = os.listdir(dir_path + dir_num)
  for file in files:
    if file.startswith('ILTJ'):
      #print (file)
      return file
    
def create_reg(list_num, dir_path, color, output_name):
  region_list = []
  str_template = 'fk5\npoint({:f},{:f})# color=' + str(color) + ' width=1 text=""'
  for num in list_num:
    file = num_filename(num, dir_path)
    [RA, DEC] = filename_to_skycoord(file)
    region_list.append((RA.deg,DEC.deg))
    #region_strs = list(map(lambda pos: 'fk5\npoint({:f},{:f})# color=red width=1 text=""'.format(*pos), region_list))
    region_strs = list(map(lambda pos: str_template.format(*pos), region_list))
  #print (region_strs)
  parser = DS9Parser('\n'.join(list(region_strs)))
  regions = parser.shapes.to_regions() 
  write_ds9(regions, output_name)
  print (output_name + " saved." )

dir_path = '/net/rijn2/data2/Haoyang/ALICE/selfcal_1b1/'
## Good selfcal, chose merged_selfcalcyle008*.h5
list1 = [8,12,16,19,22,24,26,30,32,36,43,47,50,54,61,70,73,78,85,87,88]
## Good selfcal, chose merged_selfcalcyle003*.h5
list2 = [7,14,25,31,45,46,48,57,79,81]
## Not so good selfcal, chose merged_selfcalcyle008*.h5
list3 = [4,9,23,27,29,34,35,37,49,53,56,59,66,71,76,80,82,84,86,90]
## Not so good selfcal, chose merged_selfcalcyle003*.h5
list4 = [51,83]


create_reg(list1, dir_path, 'green', 'region1.reg')
create_reg(list2, dir_path, 'yellow', 'region2.reg')
create_reg(list3, dir_path, 'red', 'region3.reg')
create_reg(list4, dir_path, 'cyan', 'region4.reg')

    
    