import numpy as np
import os
import re
import matplotlib.pylab as plt

def read_noise(filename):
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
    file_output = str(dir_num) + '_' + ms_file + '_noise_output.txt'
    file_data = str(dir_num) + '_' + ms_file + '_noise_data.txt'
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

#filename = 'selfcal_60_107992.out'
PATH = '/net/rijn2/data2/Haoyang/ALICE/selfcal_1b1/'
save_PATH = '/net/rijn2/data2/Haoyang/ALICE/dir_choice'
markers = ['o', 'v', '^', '<', '>', 's', 'p', 'P', '*', 'X', 'D']
for filename in os.listdir(PATH):
    if filename.endswith('.out'):
        dir_num = int(re.search(r'\d+', filename).group())
        print ('Processing Direction ' + str(dir_num) + ', the output file name is ' + filename)
        ms_file, number_output = read_noise(PATH + filename)
        data = split_cycle(number_output)
        cycle_num = len(data)
        print (cycle_num, 'selfcal cycles have been performed for ' + ms_file)
        print ('Plotting the noise data \n')
        plt.figure()
        for i in range(cycle_num):
            x, y = splitevenodd(data[i])
            plt.plot(x,y)
            plt.scatter(x,y, marker = markers[i],label='Cycle ' + str(i))
            plt.xticks([1,2,3,4])
        plt.title('Selfcal RMS change of ' + ms_file)
        plt.xlabel('Deconvolution Cycle Number')
        plt.ylabel('RMS (ÂµJy)')
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        plt.savefig(save_PATH + str(dir_num) + '_' + ms_file + '_noise.png', bbox_inches='tight', dpi=300)
