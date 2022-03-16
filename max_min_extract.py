#python 3
# Author: Haoyang Ye

import numpy as np
import os
import re
import matplotlib.pylab as plt
import csv
from astropy.io import fits
from pathlib import Path
import pandas as pd
import math


def max_min_val(file_name):
    print ('Opening file ', file_name)
    hdul = fits.open(file_name)
    data = hdul[0].data
    val = np.abs(data.max()/data.min())
    hdul.close()
    return val

def collect_val(dir_num, PATH, d):
    """
    dir_num: str, number of directions
    """
    d['dir_num'][int(dir_num)] = int(dir_num)
    for filename in os.listdir(PATH):
        if filename.startswith(dir_num+'_') and filename.endswith('.fits'):
            val = max_min_val(PATH + filename)
            idx = Path(filename).stem.split('_')[-1]
            d[idx][int(dir_num)] = val
            print (filename, idx, val)
    return d

def init_d(keylist,dir_sum):
    d = {}
    for i in keylist:
        d[i] = [math.nan]*(dir_sum+1)
    return d

keylist = ['dir_num', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
PATH = str(Path().absolute()).split("\/")[0] + '/fits_collection/' # Directory of current working directory
save_PATH = str(Path().absolute()).split("\/")[0] + '/data_analysis/'
dir_sum = 100

isExist = os.path.exists(save_PATH)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(save_PATH)
    print("The new directory noise_png is created!")

d = init_d(keylist, dir_sum)
for i in range(1,dir_sum+1):
    print (i)
    d = collect_val(str(i), PATH, d)

df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in d.items() ]))

df.to_hdf(save_PATH + 'max_min_fits.h5', key='selfcal_1b1', mode='w')

df.to_csv (save_PATH + r'max_min_fits.csv', index = False, header=True)
