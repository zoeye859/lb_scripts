import numpy as np
import os
import re
import matplotlib.pylab as plt
import csv

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
    #print ('The following list is saved as ' + ms_file + '_noise_output.txt \n', noise_output)
    #print ('The following list is saved as ' + ms_file + '_noise_data.txt \n',number_output)
    #file_output = save_PATH + str(dir_num) + '_' + ms_file + '_noise_output.txt'
    #file_data = save_PATH + str(dir_num) + '_' + ms_file + '_noise_data.txt'

    #with open(file_output, 'w') as f:
    #    for item in noise_output:
    #        f.write("%s\n" % item)
    #with open(file_data, 'w') as f:
    #    for item in number_output:
    #        f.write("%s\n" % item)
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
PATH = '/Users/zoeye51/Desktop/dir_choice/Output/'
save_PATH = '/Users/zoeye51/Desktop/dir_choice/png/'
markers = ['o', 'v', '^', '<', '>', 's', 'p', 'P', '*', 'X', 'D']
ratio_output = []
for filename in os.listdir(PATH):
    if filename.endswith('.out'):
        dir_num = int(re.search(r'\d+', filename).group())
        print ('Processing Direction ' + str(dir_num) + ', the output file name is ' + filename)
        ms_file, number_output = read_noise(save_PATH, PATH + filename)
        data = split_cycle(number_output)
        origin_y = data[0][1]
        cycle_num = len(data)
        print (cycle_num, 'selfcal cycles have been performed for ' + ms_file)
        print ('Plotting the noise data \n')
        plt.figure()
        min_val = []
        for i in range(cycle_num):
            x, y = splitevenodd(data[i])
            print (y)
            min_val += [min(y)]
            plt.plot(x,y)
            plt.scatter(x,y, marker = markers[i],label='Cycle ' + str(i))
            plt.xticks([1,2,3,4])
        decrease_ratio = ["{:.2%}".format((min_val[i] - origin_y)/origin_y) for i in range(cycle_num)]
        print (decrease_ratio)
        #plt.text(1, origin_y, text1)
        plt.title('Selfcal RMS change of ' + ms_file)
        plt.xlabel('Deconvolution Cycle Number')
        plt.ylabel('RMS (µJy)')
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        plt.savefig(save_PATH + str(dir_num) + '_' + ms_file + '_noise.png', bbox_inches='tight', dpi=300)

        plt.figure()
        for i in range(cycle_num):
            plt.scatter(i, float(decrease_ratio[i].strip('%')), marker = markers[i],label='Cycle ' + str(i))
        plt.plot([i for i in range(cycle_num)], [float(decrease_ratio[i].strip('%')) for i in range(cycle_num)])
        plt.xticks([0,1,2,3,4,5,6,7,8])
        plt.title('Selfcal RMS change (%)' + ms_file)
        plt.xlabel('Selfcal Cycle Number')
        plt.ylabel('Percentage (%)')
        plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
        plt.savefig(save_PATH + str(dir_num) + '_' + ms_file + '_percent.png', bbox_inches='tight', dpi=300)
        decrease_ratio.insert(0, dir_num)
        ratio_output += [decrease_ratio]

print (ratio_output)
Details = ['Dir_num', 'Cycle 0', 'Cycle 1', 'Cycle 2', 'Cycle 3', 'Cycle 4', 'Cycle 5', 'Cycle 6', 'Cycle 7', 'Cycle 8']
with open(save_PATH + 'decrease_ratio.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(Details)
    write.writerows(ratio_output)
