#python 3
# Author: Haoyang Ye

import pandas as pd
import numpy as np
from pathlib import Path

PATH = str(Path().absolute()).split("\/")[0] # Directory of current working directory
save_PATH = PATH + '/data/'

def choose_h5(num):
    if num > 5:
        return 8
    else:
        return 4

def dir_filter(max_min_fits_file, decrease_ratio_file, ):
    df1 = pd.read_hdf('max_min_fits.h5')
    df2 = pd.read_hdf('decrease_ratio.h5')
    df1 = df1.rename(columns={"0": "C0", "1": "C1", "2": "C2", "3": "C3", "4": "C4", "5": "C5", "6": "C6", "7": "C7", "8": "C8", "9": "C9"}, errors="raise")
    df1.pop("dir_num")
    val_max = df1.max(axis = 1, skipna = False).tolist()
    ind_max = df1.idxmax(axis = 1, skipna = False).tolist()
    df1['C val_max'] = val_max
    df1['C ind_max'] = ind_max
    val_min = df2.min(axis = 1, skipna = False).tolist()
    ind_min = df2.idxmin(axis = 1, skipna = False).tolist()
    df2['val_min'] = val_min
    df2['ind_min'] = ind_min
    result = pd.concat([df2, df1], axis=1, join="inner")
    result['val_min_compare'] = result['0'] - result['val_min']
    result['h5_num'] = [choose_h5(float(i)) for i in result['ind_min']]
    #### 4 conditions
    con1 = result['C9'] > result['C0'] # max/min increase from C0
    con2 = result['ind_min'] != str(0) # the biggest decrease should not be here
    con3 = result['C val_max'] > 30
    con4 = result['val_min'] < -0.1
    #con5 = result['val_min_compare'] > 0.06
    #### filter
    select = result.loc[con1 & con2 & con3 & con4]
    h5_num_filter = select['h5_num'].astype(int).tolist()
    dir_num_filter = select['dir_num'].astype(int).tolist()
    Max_min_filter = select['C val_max'].tolist()
    result.to_hdf(save_PATH + 'result.h5', key='result', mode='w')
    select.to_hdf(save_PATH + 'select.h5', key='select', mode='w')
    return h5_num_filter, dir_num_filter, Max_min_filter
