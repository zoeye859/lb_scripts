#python 3
# Author: Haoyang Ye
import numpy as np
import os
import re
import matplotlib.pylab as plt
import csv
from pathlib import Path
import pandas as pd
import argparse
import math

def init_d(keylist,dir_sum):
    d = {}
    for i in keylist:
        d[i] = [math.nan]*(dir_sum+1)
    return d

def read_noise(save_PATH, filename):
    noise_output = []
    number_output = []
    with open(filename, 'r') as outfile:
        for line in outfile:
            if line.startswith('  input MS:'):
                ms_file = line[:-1].split('/')[-1]
                #noise_output = [line[:-1].split('/')[-1]]
                #number_output = [line[:-1].split('/')[-1]]
                break
    with open(filename, 'r') as outfile:
        for line in outfile:
            if line.startswith(' == Deconvolving ('):
                temp_str = line[:-1]
                temp_num = int(re.search(r'\d+', temp_str).group())
            if line.startswith('Estimated standard deviation of background noise:'):
                noise_output += [temp_str, line[:-1]]
                number_output += [temp_num, float(re.findall(r'[\d\.\d]+', line)[0])]
    print ('The ms is: ', ms_file, '\n')
    print ('The following list is saved as ' + ms_file + '_noise_output.txt \n', noise_output)
    print ('The following list is saved as ' + ms_file + '_noise_data.txt \n',number_output)
    file_output = save_PATH + str(dir_num) + '_' + ms_file + '_noise_output.txt'
    file_data = save_PATH + str(dir_num) + '_' + ms_file + '_noise_data.txt'
    with open(file_output, 'w') as f:
        for item in noise_output:
            f.write("%s\n" % item)
    with open(file_data, 'w') as f:
        for item in number_output:
            f.write("%s\n" % item)
    return ms_file, number_output

def splitevenodd(input_list):
    return input_list[::2], input_list[1::2]

def split_cycle(input_list):
    split_output = []
    temp_idx = 0
    for i in range(1,len(input_list)):
        if input_list[i] == 1:
            split_output = split_output + [input_list[temp_idx:i]]
            temp_idx = i
    return split_output

def dir_number(filename):
    return int(re.search(r'\d+', filename).group())

parser = argparse.ArgumentParser(description='Extract background noise level for every selfcal cycle. \n \
                                Noise values and the plots of the noise level change are saved in noise_png folder. \n \
                                Make sure you are running this python script under the directory of .out files.\n')
#imaging options
parser.add_argument('--plot', help='Make plots', default=True)
args        = vars(parser.parse_args())
plot_png    = args['plot']

PATH = str(Path().absolute()).split("\/")[0] # Directory of current working directory
markers = ['o', 'v', '^', '<', '>', 's', 'p', 'P', '*', 'X', 'D']
keylist = ['dir_num', '0', '1', '2', '3', '4', '5', '6', '7', '8']
dir_sum = 100
save_PATH = PATH + '/noise_png/'
isExist = os.path.exists(save_PATH)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(save_PATH)
    print("The new directory noise_png is created!")

df_maxmin = pd.read_hdf('../max_min_fits.h5')
ratio_table = df_maxmin.T
d = init_d(keylist, dir_sum)

for filename in os.listdir(PATH):
    if filename.endswith('.out'):
        dir_num = dir_number(filename)
        print ('Processing Direction ' + str(dir_num) + ', the output file name is ' + filename)
        ms_file, number_output = read_noise(save_PATH, PATH + '/'+ filename)
        data = split_cycle(number_output)
        origin_y = data[0][1]
        cycle_num = len(data)
        print (cycle_num, 'selfcal cycles have been performed for ' + ms_file)
        min_val = []
        for i in range(cycle_num):
            x, y = splitevenodd(data[i])
            min_val += [min(y)]
        decrease_ratio = ["{:.2%}".format((min_val[i] - origin_y)/origin_y) for i in range(cycle_num)]
        decrease_ratio_val = [(min_val[i] - origin_y)/origin_y for i in range(cycle_num)]
        print ('Decreased ratio are: ', decrease_ratio)
        d['dir_num'][int(dir_num)] = int(dir_num)
        for i in range(len(decrease_ratio_val)):
            d[str(i)][int(dir_num)] = decrease_ratio_val[i]
        df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in d.items() ]))
        df.to_hdf(save_PATH + 'decrease_ratio.h5', key='decrease_ratio', mode='w')
        df.to_csv (save_PATH + 'decrease_ratio.csv', index = False, header=True)
        if plot_png:
            print ('Plotting the noise data \n')
            plt.figure()
            min_val = []
            for i in range(cycle_num):
                x, y = splitevenodd(data[i])
                plt.plot(x, y)
                plt.scatter(x,y, marker = markers[i],label='Cycle ' + str(i))
            plt.title(str(dir_num) + '_' + ms_file)
            plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
            plt.savefig(save_PATH + str(dir_num) + '_' + ms_file + '_noise.png', bbox_inches='tight', dpi=300)

            plt.figure()
            fig, ax1 = plt.subplots()
            color = 'tab:red'
            ax1.set_xlabel('Selfcal Cycle Number', color=color)
            ax1.set_ylabel('Noise Decrease Percentage (%)', color=color)
            ax1.tick_params(axis='y', labelcolor=color)
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            color1 = 'tab:blue'
            ax2.set_ylabel('Max/Min', color=color1)  # we already handled the x-label with ax1
            ax2.tick_params(axis='y', labelcolor=color1)
            for i in range(cycle_num):
                ax1.scatter(i, float(decrease_ratio[i].strip('%')), marker = markers[i],label='Cycle ' + str(i))
                ax2.scatter(i, ratio_table[int(dir_num)][i], marker = markers[i],label='Cycle ' + str(i))
            ax1.plot([i for i in range(cycle_num)], [float(decrease_ratio[i].strip('%')) for i in range(cycle_num)], c='r')
            ax2.plot([i for i in range(cycle_num)], [ratio_table[int(dir_num)][i] for i in range(cycle_num)], c='b')
            plt.title(str(dir_num) + '_' + ms_file)
            ax1.legend(bbox_to_anchor=(1.08,1), loc="upper left")
            plt.savefig(save_PATH + str(dir_num) + '_' + ms_file + '_percent.png', bbox_inches='tight', dpi=300)
