import os

dir_path = '/net/rijn2/data2/Haoyang/ALICE/selfcal_1b1/'
os.chdir(dir_path)
lines = []

#90 directories
for i in range(1,91):
  print (i)
  lines += [str(i) + '\n']
  file_path = dir_path + str(i)
  files = os.listdir(file_path)
  for file in files:
    if file.startswith('ILTJ'):
      print (file)
      lines += [file + '\n']

#print (lines)
os.chdir('/net/rijn2/data2/Haoyang/ALICE/1arcsec')
with open('dir_name_list.txt', 'w') as f:
  for line in lines:
    f.write(line)
    
