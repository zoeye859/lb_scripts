# !/usr/bin/python
# -*- coding: utf-8 -*-

import os
import numpy as np
import tables
import shutil

# The selfcal directions are in individual folders now

def collect_png_bynumber(dir_path, cycle_num=10, png_dir_name = 'png_collection'):
    '''
    After running selfcal on multiple directions chosen from 6'' image, we gathered multiple folders.
    Each folder contains the selfcal results for one direction.
    We want to collect all the .png files from each of these folders to exam.
    Args:
        dir_path (str): the path to the directory in which the selfcal results are.
        png_dir_name (str): name of the png collection directory, default: png_collection.
        cycle_num (int): maximum cycle number of the self-cal, default: 10
    Returns:
        png_output_path (str): the path to the png collection directory where we put all the .png files
    '''
    dir_path = str(dir_path)
    dir_name = [d for d in os.listdir(dir_path) if os.path.isdir(d)]
    png_dir_name = str(png_dir_name)
    os.chdir(dir_path)
    os.mkdir(png_dir_name)
    png_output_path = os.path.join(dir_path, png_dir_name)
    for name in dir_name:
        os.chdir(os.path.join(dir_path, name)) # Enter each direction folder
        print ('Copying png files from directory ', os.getcwd())
        for i in range(cycle_num):
            png_name = name + '_' + str(i) + '.png'
            if os.path.isfile('image_'+str(i)+'.png'): 
                print ('Copying file name:', png_name)
                shutil.copy('image_'+str(i)+'.png', png_output_path+'/'+png_name) 
    print ('The number of collected png file:', ?????)
    return png_output_path

def collect_fits_bynumber(dir_path, cycle_num=10, fits_dir_name = 'fits_collection'):
    '''
    After running selfcal on multiple directions chosen from 6'' image, we gathered multiple folders.
    Each folder contains the selfcal results for one direction.
    We want to collect all the .fits files from each of these folders to exam.
    Args:
        dir_path (str): the path to the directory in which the selfcal results are.
        fits_dir_name (str): Name of the fits collection directory, default: fits_collection.
        cycle_num (int): maximum cycle number of the self-cal, default: 10
    Returns:
        fits_output_path (str): the path to the fits collection directory where we put all the .fits files
    '''
    dir_path = str(dir_path)
    dir_name = [d for d in os.listdir(dir_path) if os.path.isdir(d)]
    fits_dir_name = str(fits_dir_name)
    os.chdir(dir_path)
    os.mkdir(fits_dir_name)
    fits_output_path = os.path.join(dir_path, fits_dir_name)
    for name in dir_name:
        os.chdir(os.path.join(dir_path, name)) # Enter each direction folder
        print ('Copying fits files from directory ', os.getcwd())
        for i in range(cycle_num):
            fits_name = name + '_' + str(i) + '.fits'
            if os.path.isfile('image_00'+str(i)+'-MFS-image.fits'): 
                print ('Copying file name:', fits_name)
                shutil.copy('image_00'+str(i)+'-MFS-image.fits', fits_output_path+'/'+fits_name) 
    return fits_output_path

# make a directory to gather all png files

# choose .h5 files

# random test

# cs stations

# test

# copy sol001 to sol000, delete sol001, or use Jurjen's option

# h5 merger

# make screen


def main(dir_path, png_dir_name = 'png_collection'):
    '''
    
    
    
    Args:
        dir_path (str): the path to the directory in which the selfcal results are.
        png_dir_name (str): Name of the png collection directory, default: png_collection.
    Returns:
        
    '''
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='From selfcal results to a screen for WSCLEAN and IDG to use')

    parser.add_argument('h5parm', type=str,
                        help='H5parm to which this action should be performed .')
    parser.add_argument('MSfiles', type=str,
                        help='MS for which the new solset shall be created.')
    parser.add_argument('--solset_in', type=str, default='sol000',
                        help='Input solution set')
    parser.add_argument('--solset_out', type=str, default='sol001',
                        help='Output solution set (has to be different from input solution set)')
    parser.add_argument('--soltab_list', '--soltab_list', type=str, default='phase000,amplitude000',
                        help='Comma-separated list of soltabs to be copied')
    parser.add_argument('--superstation', type=str, default='ST001',
                        help='Reference station from which data should be copied')
    parser.add_argument('--restrictToCS',
                        help='Restrict the copy action to core stations only',
                        action='store_true', dest="restrictToCS")
    parser.add_argument('--match_ptg',
			help='Match the pointing direction of the MS, for transferring between directions',
			action='store_true', dest='match_ptg')
    parser.add_argument('--useh5coordinates', help='Copy the source coordinate from the input .h5 file at sol000', default=True)

    args = parser.parse_args()

    format_stream = logging.Formatter("%(asctime)s\033[1m %(levelname)s:\033[0m %(message)s", "%Y-%m-%d %H:%M:%S")
    format_file = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    logging.root.setLevel(logging.INFO)

    log = logging.StreamHandler()
    log.setFormatter(format_stream)
    logging.root.addHandler(log)

    MSfiles = args.MSfiles
    h5parmfile = args.h5parm

    main(h5parmfile, MSfiles, solset_in=args.solset_in, solset_out=args.solset_out, soltab_list=args.soltab_list, superstation=args.superstation, restrictToCS=args.restrictToCS, matchPtg=args.match_ptg, useh5coordinates=args.useh5coordinates)
