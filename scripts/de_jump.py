# !/usr/bin/python
# -*- coding: utf-8 -*-

import os
import numpy as np
import shutil

dir_path = '/net/rijn2/data2/Haoyang/ALICE/selfcal_1b1/'
parset_path = '/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_2022jan/jump_check/losoto_plottec.parset'
paste_path1='/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_2022jan/jump_check/green'
paste_path2='/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_2022jan/jump_check/yellow'
paste_path3 = '/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_2022jan/jump_check/puzzle'
paste_path4 = '/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_2022jan/jump_check/red'
paste_path5 = '/net/rijn2/data2/Haoyang/ALICE/selected_h5/h5_collection_2022jan/jump_check/cyan'

for i in [8,12,16,19,22,24,26,30,32,36,50,54,61,70,73,78,85,87,88]:
  print ("Copying h5 file from folder " + str(i))
  os.chdir(dir_path+str(i))
  for file in os.listdir("."):
    if file.startswith('tec0_selfcalcyle008') and file.endswith('_concat.ms.avg.h5'):
        shutil.copy(file, paste_path1)
        file_name = str(i)+'_tec0_selfcalcyle008.h5'
        file_name2 = file
        print ("Copying " + file_name)
        shutil.copy(file, paste_path1+'/'+ file_name) 
  os.chdir(paste_path1)
  cmd = 'losoto '+file_name+ ' ' + parset_path 
  os.system(cmd)
  #shutil.copy('plottec/tecdirpointing_freq144627380.37109375.png', 'plottec/'+str(i)+'tecdirpointing.png')
  shutil.copy('plottec/tecdirpointing_freq144627380.37109375.png', 'plottec/'+file_name2+'.png')



for i in [7,25,31,45,46,48,57,79,81]:
  print ("Copying h5 file from folder " + str(i))
  os.chdir(dir_path+str(i))
  for file in os.listdir("."):
    if file.startswith('tec0_selfcalcyle003') and file.endswith('_concat.ms.avg.h5'):
        shutil.copy(file, paste_path2)
        file_name = str(i)+'_tec0_selfcalcyle003.h5'
        file_name2 = file
        print ("Copying " + file_name)
        shutil.copy(file, paste_path2+'/'+ file_name)
  os.chdir(paste_path2)
  cmd = 'losoto '+file_name+ ' ' + parset_path 
  os.system(cmd)
  shutil.copy('plottec/tecdirpointing_freq144627380.37109375.png', 'plottec/'+file_name2+'.png')

for i in [43, 47]:
  print ("Copying h5 file from folder " + str(i))
  os.chdir(dir_path+str(i))
  for file in os.listdir("."):
    if file.startswith('tec0_selfcalcyle008') and file.endswith('_concat.ms.avg.h5'):
        shutil.copy(file, paste_path3)
        file_name = str(i)+'_tec0_selfcalcyle008.h5'
        file_name2 = file
        print ("Copying " + file_name)
        shutil.copy(file, paste_path3+'/'+ file_name)
  os.chdir(paste_path3)
  cmd = 'losoto '+file_name+ ' ' + parset_path 
  os.system(cmd)
  shutil.copy('plottec/tecdirpointing_freq144627380.37109375.png', 'plottec/'+file_name2+'.png')

for i in [14]:
  print ("Copying h5 file from folder " + str(i))
  os.chdir(dir_path+str(i))
  for file in os.listdir("."):
    if file.startswith('tec0_selfcalcyle003') and file.endswith('_concat.ms.avg.h5'):
        shutil.copy(file, paste_path3)
        file_name = str(i)+'_tec0_selfcalcyle003.h5'
        file_name2 = file
        print ("Copying " + file_name)
        shutil.copy(file, paste_path3+'/'+ file_name)
  os.chdir(paste_path3)
  cmd = 'losoto '+file_name+ ' ' + parset_path 
  os.system(cmd)
  shutil.copy('plottec/tecdirpointing_freq144627380.37109375.png', 'plottec/'+file_name2+'.png')

for i in [4,9,23,27,29,34,35,37,49,53,56,59,66,71,76,80,82,84,86,90]:
  print ("Copying h5 file from folder " + str(i))
  os.chdir(dir_path+str(i))
  for file in os.listdir("."):
    if file.startswith('tec0_selfcalcyle008') and file.endswith('_concat.ms.avg.h5'):
        shutil.copy(file, paste_path4)
        file_name = str(i)+'_tec0_selfcalcyle008.h5'
        file_name2 = file
        print ("Copying " + file_name)
        shutil.copy(file, paste_path4+'/'+ file_name)
  os.chdir(paste_path4)
  cmd = 'losoto '+file_name+ ' ' + parset_path 
  os.system(cmd)
  shutil.copy('plottec/tecdirpointing_freq144627380.37109375.png', 'plottec/'+file_name2+'.png')

for i in [51, 83]:
  print ("Copying h5 file from folder " + str(i))
  os.chdir(dir_path+str(i))
  for file in os.listdir("."):
    if file.startswith('tec0_selfcalcyle003') and file.endswith('_concat.ms.avg.h5'):
        shutil.copy(file, paste_path5)
        file_name = str(i)+'_tec0_selfcalcyle003.h5'
        file_name2 = file
        print ("Copying " + file_name)
        shutil.copy(file, paste_path5+'/'+ file_name)
  os.chdir(paste_path5)
  cmd = 'losoto '+file_name+ ' ' + parset_path 
  os.system(cmd)
  shutil.copy('plottec/tecdirpointing_freq144627380.37109375.png', 'plottec/'+file_name2+'.png')
