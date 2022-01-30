#python 3
# Author: Haoyang Ye

from astropy.io import fits
import tables
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from matplotlib.patches import Rectangle
import re
import os
from pathlib import Path
from shutil import copyfile, rmtree

PATH = str(Path().absolute()).split("\/")[0] # Directory of current working directory
dir_path = PATH + '/screen_plots/'
ext_scale = 1.5 #make the screen 1.5 times bigger than the given boxwidth

def plot_screen_wcs(fits_data, fits_header, size, antenna, time_slot = 10, freq_slot = 2):
    '''
    plot the screen .fits image with wcs axises

    fits_data: numpy array taken from .fits file
    fits_header: the header of the .fits file, to take the phase centre values
    boxwidth: float, the size of the screen in both axises in the unit of degrees
    antenna: str list, output from function get_antenna_names
    time_slot: int, default 11th time_slot
    freq_slot: int, default 3th frequency slot
    '''
    if os.path.exists(dir_path)==False:
        os.mkdir(dir_path)
    fig = plt.figure(figsize=(21, 9))
    wcs = WCS(fits_header, naxis = 2)
    image_data = fits_data[time_slot,freq_slot,:,:,:,:]
    title_list = ['XX real', 'XX imaginary', 'YY real', 'YY imaginary', 'XX phase',
     'XX amplitude', 'YY phase', 'YY amplitude']
    position = (size/2, size/2) # image centre
    cutout = Cutout2D(image_data[0,0,:,:], position, (size/ext_scale, size/ext_scale))# show the boxwidth size
    dir_x_pix, dir_y_pix, dir_tag_num = deg2pix(h5file, boxwidth, size, centre_RA, centre_DEC) # convert direction into pixels

    ## plot the screen for one of the core stations, all core stations' screen should be
    #identical
    ## Here I plot the first station
    print ("Plotting the screen for core stations with wcs axises")

    ## plot XX real, XX imag, YY real, YY imag
    for i in range(4):
        ax1 = plt.subplot(2,4,i+1,projection=wcs)
        ax1.set_title(title_list[i])
        im1 = ax1.imshow(image_data[0,i,:,:], cmap=plt.cm.viridis,
             interpolation='nearest', origin='lower')
        ax1.set_xlabel('Right Ascension (J2000)')
        ax1.set_ylabel('Declination (J2000)')
        ax1.grid(alpha=0.6)
        ax1.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
        for idx, txt in enumerate(dir_tag_num):
            ax1.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="r")
        fig.colorbar(im1)
        #show a square of the given boxwidth
        cutout.plot_on_original(color='white')

    ## plot 'XX phase', 'XX amplitude', 'YY phase', 'YY amplitude'
    for i in [0, 2]:
        ax2 = plt.subplot(2,4,i+6,projection=wcs)
        ax2.set_title(title_list[i+5])
        im2 = ax2.imshow(np.sqrt(image_data[0,i,:,:]**2 + image_data[0,i+1,:,:]**2), cmap=plt.cm.viridis,
             interpolation='nearest', origin='lower')
        ax2.set_xlabel('Right Ascension (J2000)')
        ax2.set_ylabel('Declination (J2000)')
        ax2.grid(alpha=0.6)
        ax2.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
        for idx, txt in enumerate(dir_tag_num):
            ax2.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="r")
        fig.colorbar(im2)
        cutout.plot_on_original(color='white')

        ax3 = plt.subplot(2,4,i+5,projection=wcs)
        ax3.set_title(title_list[i+4])
        im3 = ax3.imshow(np.angle(image_data[0,i,:,:] + image_data[0,i+1,:,:]*1j), cmap=plt.cm.viridis,
             interpolation='nearest', origin='lower')
        ax3.set_xlabel('Right Ascension (J2000)')
        ax3.set_ylabel('Declination (J2000)')
        ax3.grid(alpha=0.6)
        ax3.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
        for idx, txt in enumerate(dir_tag_num):
            ax3.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="r")
        fig.colorbar(im3)
        cutout.plot_on_original(color='white')

    fig.suptitle('Core Station')
    plt.savefig(dir_path + 'core_station_screen_wcs.png', bbox_inches='tight')
    #plt.show()
    plt.close()

    for ant_idx in range(len(antenna)-1):
        fig = plt.figure(figsize=(21, 9))
        ant_name = antenna[-2-ant_idx]
        print ("Plotting the screen for " + ant_name)
        ## plot the screen for each of the remote/international stations
        data_idx = -1-ant_idx

        ## plot XX real, XX imag, YY real, YY imag
        for i in range(4):
            ax2 = plt.subplot(2,4,i+1,projection=wcs)
            ax2.set_title(title_list[i])
            im2 = ax2.imshow(image_data[data_idx,i,:,:], cmap='gnuplot2',
                 interpolation='nearest', origin='lower')
            ax2.set_xlabel('Right Ascension (J2000)')
            ax2.set_ylabel('Declination (J2000)')
            ax2.grid(alpha=0.6)
            ax2.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
            for idx, txt in enumerate(dir_tag_num):
                ax2.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="cyan")
            fig.colorbar(im2)
            cutout.plot_on_original(color='white')

        ## plot 'XX phase', 'XX amplitude', 'YY phase', 'YY amplitude'
        for i in [0, 2]:
            ax3 = plt.subplot(2,4,i+6,projection=wcs)
            ax3.set_title(title_list[i+5])
            im3 = ax3.imshow(np.sqrt(image_data[data_idx,i,:,:]**2 + image_data[data_idx,i+1,:,:]**2), cmap='viridis',
                 interpolation='nearest', origin='lower')
            #print (np.sqrt(image_data[data_idx,i,:,:]**2 + image_data[data_idx,i+1,:,:]**2))
            ax3.set_xlabel('Right Ascension (J2000)')
            ax3.set_ylabel('Declination (J2000)')
            ax3.grid(alpha=0.6)
            ax3.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
            for idx, txt in enumerate(dir_tag_num):
                ax3.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="cyan")
            fig.colorbar(im3)
            cutout.plot_on_original(color='white')

            ax4 = plt.subplot(2,4,i+5,projection=wcs)
            ax4.set_title(title_list[i+4])
            im4 = ax4.imshow(np.angle(image_data[data_idx,i,:,:] + image_data[data_idx,i+1,:,:]*1j), cmap='hsv',
                 interpolation='nearest', origin='lower')
            ax4.set_xlabel('Right Ascension (J2000)')
            ax4.set_ylabel('Declination (J2000)')
            ax4.grid(alpha=0.6)
            ax4.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
            for idx, txt in enumerate(dir_tag_num):
                ax4.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="cyan")
            fig.colorbar(im4)
            cutout.plot_on_original(color='white')

        fig.suptitle(ant_name)
        plt.savefig(dir_path + ant_name + '_screen_wcs.png', bbox_inches='tight')
        #plt.show()
        plt.close()



def plot_screen(fits_data, fits_header, size, antenna, time_slot = 10, freq_slot = 2):
    '''
    plot the screen .fits image without wcs axises

    fits_data: numpy array taken from .fits file
    fits_header: the header of the .fits file, to take the phase centre values
    boxwidth: float, the size of the screen in both axises in the unit of degrees
    antenna: str list, output from function get_antenna_names
    time_slot: int, default 11th time_slot
    freq_slot: int, default 3th frequency slot
    '''
    if os.path.exists(dir_path)==False:
        os.mkdir(dir_path)
    fig = plt.figure(figsize=(21, 9))
    image_data = fits_data[time_slot,freq_slot,:,:,:,:]
    title_list = ['XX real', 'XX imaginary', 'YY real', 'YY imaginary', 'XX phase',
     'XX amplitude', 'YY phase', 'YY amplitude']
    dir_x_pix, dir_y_pix, dir_tag_num = deg2pix(h5file, boxwidth, size, centre_RA, centre_DEC) # convert direction into pixels

    ## plot the screen for one of the core stations, all core stations' screen should be
    #identical
    ## Here I plot the first station
    print ("Plotting the screen for core stations with pixel axises")

    ## plot XX real, XX imag, YY real, YY imag
    for i in range(4):
        ax1 = plt.subplot(2,4,i+1)
        ax1.set_title(title_list[i])
        im1 = ax1.imshow(image_data[0,i,:,:], cmap=plt.cm.viridis,
             interpolation='nearest', origin='lower')
        ax1.set_xlabel('RA (pixel))')
        ax1.set_ylabel('DEC (pixel)')
        ax1.grid(alpha=0.6)
        ax1.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
        for idx, txt in enumerate(dir_tag_num):
            ax1.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="r")
        fig.colorbar(im1)
        #show a square of the given boxwidth
        rect = Rectangle((size//2 - size/ext_scale/2,size//2 - size/ext_scale/2),size/ext_scale,size/ext_scale,linewidth=1,edgecolor='w',facecolor='none', alpha = 0.8)
        ax1.add_patch(rect)
        #ax1.invert_xaxis()


    ## plot 'XX phase', 'XX amplitude', 'YY phase', 'YY amplitude'
    for i in [0, 2]:
        ax2 = plt.subplot(2,4,i+6)
        ax2.set_title(title_list[i+5])
        im2 = ax2.imshow(np.sqrt(image_data[0,i,:,:]**2 + image_data[0,i+1,:,:]**2), cmap=plt.cm.viridis,
             interpolation='nearest', origin='lower')
        ax2.set_xlabel('RA (pixel))')
        ax2.set_ylabel('DEC (pixel)')
        ax2.grid(alpha=0.6)
        ax2.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
        for idx, txt in enumerate(dir_tag_num):
            ax2.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="r")
        fig.colorbar(im2)
        rect = Rectangle((size//2 - size/ext_scale/2,size//2 - size/ext_scale/2),size/ext_scale,size/ext_scale,linewidth=1,edgecolor='w',facecolor='none', alpha = 0.8)
        ax2.add_patch(rect)
        #ax2.invert_xaxis()


        ax3 = plt.subplot(2,4,i+5)
        ax3.set_title(title_list[i+4])
        im3 = ax3.imshow(np.angle(image_data[0,i,:,:] + image_data[0,i+1,:,:]*1j), cmap=plt.cm.viridis,
             interpolation='nearest', origin='lower')
        ax3.set_xlabel('RA (pixel))')
        ax3.set_ylabel('DEC (pixel)')
        ax3.grid(alpha=0.6)
        ax3.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
        for idx, txt in enumerate(dir_tag_num):
            ax3.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="r")
        fig.colorbar(im3)
        rect = Rectangle((size//2 - size/ext_scale/2,size//2 - size/ext_scale/2),size/ext_scale,size/ext_scale,linewidth=1,edgecolor='w',facecolor='none', alpha = 0.8)
        ax3.add_patch(rect)
        #ax3.invert_xaxis()

    fig.suptitle('Core Station')
    plt.savefig(dir_path + 'core_station_screen.png', bbox_inches='tight')
    #plt.show()
    plt.close()

    for ant_idx in range(len(antenna)-1):
        fig = plt.figure(figsize=(21, 9))
        ant_name = antenna[-2-ant_idx]
        print ("Plotting the screen for " + ant_name)
        ## plot the screen for each of the remote/international stations
        data_idx = -1-ant_idx

        ## plot XX real, XX imag, YY real, YY imag
        for i in range(4):
            ax2 = plt.subplot(2,4,i+1)
            ax2.set_title(title_list[i])
            im2 = ax2.imshow(image_data[data_idx,i,:,:], cmap='gnuplot2',
                 interpolation='nearest', origin='lower')
            ax2.set_xlabel('RA (pixel))')
            ax2.set_ylabel('DEC (pixel)')
            ax2.grid(alpha=0.6)
            ax2.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
            for idx, txt in enumerate(dir_tag_num):
                ax2.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="cyan")
            fig.colorbar(im2)
            rect = Rectangle((size//2 - size/ext_scale/2,size//2 - size/ext_scale/2),size/ext_scale,size/ext_scale,linewidth=1,edgecolor='w',facecolor='none', alpha = 0.8)
            ax2.add_patch(rect)
            #ax2.invert_xaxis()

        ## plot 'XX phase', 'XX amplitude', 'YY phase', 'YY amplitude'
        for i in [0, 2]:
            ax3 = plt.subplot(2,4,i+6)
            ax3.set_title(title_list[i+5])
            im3 = ax3.imshow(np.sqrt(image_data[data_idx,i,:,:]**2 + image_data[data_idx,i+1,:,:]**2), cmap=plt.cm.viridis,
                 interpolation='nearest', origin='lower')
            #print (np.sqrt(image_data[data_idx,i,:,:]**2 + image_data[data_idx,i+1,:,:]**2))
            ax3.set_xlabel('RA (pixel))')
            ax3.set_ylabel('DEC (pixel)')
            ax3.grid(alpha=0.6)
            ax3.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
            for idx, txt in enumerate(dir_tag_num):
                ax3.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="cyan")
            fig.colorbar(im3)
            rect = Rectangle((size//2 - size/ext_scale/2,size//2 - size/ext_scale/2),size/ext_scale,size/ext_scale,linewidth=1,edgecolor='w',facecolor='none', alpha = 0.8)
            ax3.add_patch(rect)
            #ax3.invert_xaxis()

            ax4 = plt.subplot(2,4,i+5)
            ax4.set_title(title_list[i+4])
            im4 = ax4.imshow(np.angle(image_data[data_idx,i,:,:] + image_data[data_idx,i+1,:,:]*1j), cmap='hsv',
                 interpolation='nearest', origin='lower')
            ax4.set_xlabel('RA (pixel))')
            ax4.set_ylabel('DEC (pixel)')
            ax4.grid(alpha=0.6)
            ax4.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 5, edgecolor = 'r')
            for idx, txt in enumerate(dir_tag_num):
                ax4.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="cyan")
            fig.colorbar(im4)
            rect = Rectangle((size//2 - size/ext_scale/2,size//2 - size/ext_scale/2),size/ext_scale,size/ext_scale,linewidth=1,edgecolor='w',facecolor='none', alpha = 0.8)
            ax4.add_patch(rect)
            #ax4.invert_xaxis()

        fig.suptitle(ant_name)
        plt.savefig(dir_path + ant_name + '_screen.png', bbox_inches='tight')
        #plt.show()
        plt.close()


def plot_directions(h5file, centre_RA, centre_DEC, boxwidth):
    '''
    plot the directions within the wcs framework

    fits_header: the header of the .fits file, to take the phase centre values
    h5file: str, input merged h5parm name
    direction_array: numpy array, output of the get_direction function
    '''
    dir_array, dir_tag, dir_x, dir_y, dir_x_deg, dir_y_deg = get_direction(h5file)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(dir_x_deg, dir_y_deg)
    for i, txt in enumerate(dir_tag):
        plt.annotate(txt, (dir_x_deg[i], dir_y_deg[i]))
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title('Direction Plot')
    ax.scatter(centre_RA, centre_DEC)
    plt.annotate('centre', (centre_RA, centre_DEC))
    #draw_circle1 = plt.Circle((centre_RA, centre_DEC), boxwidth/2, fill=False)
    #draw_circle2 = plt.Circle((centre_RA, centre_DEC), boxwidth/2*ext_scale, fill=False)
    ax.set_xlim([centre_RA - boxwidth, centre_RA + boxwidth])
    ax.set_ylim([centre_DEC - boxwidth, centre_DEC + boxwidth])
    ax.set_aspect(1)
    #ax.add_artist(draw_circle1)
    #ax.add_artist(draw_circle2)
    rect1 = Rectangle((centre_RA - boxwidth*ext_scale/2,centre_DEC - boxwidth*ext_scale/2),boxwidth*ext_scale,boxwidth*ext_scale,linewidth=1,edgecolor='b',facecolor='none', alpha = 0.5)
    rect2 = Rectangle((centre_RA - boxwidth/2,centre_DEC - boxwidth/2),boxwidth,boxwidth,linewidth=1,edgecolor='b',facecolor='none', alpha = 0.2)
    ax.add_patch(rect1)
    ax.add_patch(rect2)
    ax.invert_xaxis()
    #ax.annotate(str(boxwidth)+' region', (centre_RA,centre_DEC + boxwidth/2))
    #ax.annotate(str(boxwidth)+'1.5 region', (centre_RA,centre_DEC + boxwidth*1.5/2))
    plt.savefig(dir_path + 'Direction_plot.png', bbox_inches='tight')
    plt.show()


def select_direction2remove(dir_num, remove_num, random_t = True, remove_list_manual = []):
    '''
    select directions to remove randomly or manually

    dir_num: int, direction number taken from the .h5 file
    remove_num: int, the number of directions want to be removed
    random_t: bool, true for random, false for manual
    remove_list_manual: int list, if want to remove directions manually, provide a list of direction numbers
    '''
    if remove_num > dir_num:
        raise Exception("You are removing more directions than existing directions")
    if random_t == True:
        remove_list = random.sample(range(0, dir_num), remove_num)
    elif random_t == False:
        remove_list = remove_list_manual
    print ('You are removing ' + str(len(remove_list)) + ' directions, they are: ', remove_list)
    return remove_list

def deg2pix(h5file, boxwidth, size, centre_RA, centre_DEC):
    '''
    convert degrees to pixels in the screen image plane
    dir_x_pix: float list - dirction RAs in pixels
    dir_y_pix: float list - dirction DECs in pixels
    dir_tag_num: string list - direction tags only in number
    '''
    dir_array, dir_tag, dir_x, dir_y, dir_x_deg, dir_y_deg = get_direction(h5file)
    ref_pos_RA = centre_RA + boxwidth * ext_scale/2
    ref_pos_DEC = centre_DEC - boxwidth * ext_scale/2
    dir_x_pix = [(ref_pos_RA-x)* size/ext_scale/boxwidth for x in dir_x_deg]
    dir_y_pix = [(y - ref_pos_DEC)* size/ext_scale/boxwidth for y in dir_y_deg]
    dir_tag_num = [re.split('(\d+)', tag)[1] for tag in dir_tag]
    return dir_x_pix, dir_y_pix, dir_tag_num

def get_antenna_names(h5file):
    '''
    plot directions over the screen .fits image

    h5file: str, input merged h5parm name
    antenna: string list, the name of all antennas read from the h5 parm
    e.g.
    ['RS106HBA', 'RS205HBA', 'RS208HBA', 'RS210HBA', 'RS305HBA', 'RS306HBA', 'RS307HBA',
    'RS310HBA', 'RS406HBA', 'RS407HBA', 'RS409HBA', 'RS503HBA', 'RS508HBA', 'RS509HBA',
    'DE601HBA', 'DE602HBA', 'DE603HBA', 'DE604HBA', 'DE605HBA', 'FR606HBA', 'SE607HBA',
    'UK608HBA', 'DE609HBA', 'PL610HBA', 'PL611HBA', 'PL612HBA', 'IE613HBA', 'ST001']

    '''
    H = tables.open_file(h5file)
    antenna_list = H.root.sol000.antenna[:]
    antenna_num = len(antenna_list)
    print ('There are ' + str(antenna_num) + 'remote/international stations.\n')
    print ('The input h5parm contains the following antenna info:\n', antenna_list)
    antenna = [antenna_list[i][0] for i in range(antenna_num)]
    H.close()
    return [str(ant).split("\'")[1] for ant in antenna]

def get_field_centre(fits_header):
    '''
    obtain the field centre from the .fits image

    fits_header: str, input merged h5parm name
    centre_RA: float, RA in decimal degrees
    centre_DEC: float, DEC in decimal degrees
    e.g.
    (203.588312, 37.323366)
    '''
    centre_RA = fits_header['CRVAL1'] # RA of reference point, unit: deg
    centre_DEC = fits_header['CRVAL2'] # DEC of reference point, unit: deg
    return centre_RA, centre_DEC

def get_direction(h5file):
    '''
    obtain the directions from the .h5 file

    h5file: str, input merged h5parm name
    direction_array: numpy array, name and (RA, DEC) in radians of the directions
    e.g
    array([('Dir00', [-2.7129,  0.6202]), ('Dir01', [-2.7331,  0.6322]),
       ('Dir02', [-2.7485,  0.6765]), ('Dir03', [-2.7078,  0.6676]),
       ('Dir04', [-2.6872,  0.6464]), ('Dir05', [-2.6935,  0.642 ]),
       ('Dir06', [-2.748 ,  0.6248]), ('Dir07', [-2.6978,  0.6422]),
       ('Dir08', [-2.7042,  0.6438]), ('Dir09', [-2.7437,  0.6734]),
       ('Dir10', [-2.7234,  0.6807]), ('Dir11', [-2.7026,  0.6285]),
       ('Dir12', [-2.7107,  0.644 ])],
      dtype=[('name', 'S128'), ('dir', '<f4', (2,))])

    '''
    H = tables.open_file(h5file)
    dir_array = H.root.sol000.source[:]
    dir_num = len(dir_array)
    #print ('There are ' + str(dir_num) + ' directions.')
    #print ('These directions are:\n', dir_array)
    H.close()
    #print (dir_array[1][1][0], dir_array[1][1][1], dir_array[1][0])
    dir_tag = [str(dir_array[i][0]).split("\'")[1] for i in range(dir_num)]
    dir_x = [dir_array[idx][1][0] for idx in range(dir_num)]
    dir_y = [dir_array[idy][1][1] for idy in range(dir_num)]
    dir_x_deg = [np.rad2deg(x) if x>0 else np.rad2deg(x+2*np.pi) for x in dir_x]
    dir_y_deg = [np.rad2deg(y) for y in dir_y]
    dir_xy = np.column_stack((dir_tag, dir_x_deg, dir_y_deg))
    #print ('In degrees are:\n', dir_xy)
    return dir_array, dir_tag, dir_x, dir_y, dir_x_deg, dir_y_deg

def dir_h5_file(h5file_ind):
    '''
    obtain the directions from the .h5 file

    h5file: str, input merged h5parm name
    h5file_ind: the starting letters of the individual h5 files,
    to identify individual h5 files from the merged h5 file example 'merged_selfcalcyle'
    '''
    dir_h5 = []
    print ('Finding the individual .h5 files in the current directory...')
    for filename in os.listdir(PATH):
        if filename.startswith(str(h5file_ind)) and filename.endswith('.ms.avg.h5'):
            H = tables.open_file(filename)
            dir_now = H.root.sol000.source[:][0][1]
            dir_num = np.where(dir_xy_pos == dir_now)[0][0]
            ra = deg2HMS(dir_x_deg[dir_num])
            dec = deg2DMS(dir_y_deg[dir_num])
            dir_h5 += [[dir_num, ra + '+' + dec, filename]]
            H.close()
    textfile = open(dir_path + "h5_names.txt", "w")
    for element in dir_h5:
        textfile.write(str(element) + "\n")
    textfile.close()
    return dir_h5

def deg2DMS(deg_input, decimal_num=2):
    '''
    decimal degree to degree arcminute arcsecond

    deg_input: float, input decimal degrees
    decimal_num: int, number of decimals, default 2
    '''
    deg = np.floor(deg_input)
    arc_minute = np.floor(60*(deg_input - deg))
    arc_second = (60*(deg_input - deg) - arc_minute)*60
    arc_second_int = np.floor(arc_second)
    arc_second_dec = arc_second - arc_second_int
    if decimal_num == 0:
        return '%02d%02d%02d' % (deg, arc_minute, np.round(arc_second))
    elif decimal_num < 0:
        raise Exception("The given decimal number has to be equal or larger than 0")
    else:
        format_str = '%0' + str(decimal_num) + 'd'
        arc_second_dec *= 10**decimal_num
        return '%02d%02d%02d' % (deg, arc_minute, arc_second_int) + '.' + format_str % (arc_second_dec)


def deg2HMS(deg_input, decimal_num=2):
    '''
    decimal degree to hour minute second

    deg_input: float, input decimal degrees
    decimal_num: int, number of decimals, default 2
    '''
    hour = int(deg_input/15)
    minute = np.floor((deg_input/15-hour)*60)
    second = ((deg_input/15-hour)*60-minute)*60
    second_int = np.floor(second)
    second_dec = second - second_int
    if decimal_num == 0:
        return '%02d%02d%02d' % (hour, minute, np.round(second))
    elif decimal_num < 0:
        raise Exception("The given decimal number has to be equal or larger than 0")
    else:
        format_str = '%0' + str(decimal_num) + 'd'
        second_dec *= 10**decimal_num
        return '%02d%02d%02d' % (hour, minute, second_int) + '.' + format_str % (second_dec)


parser = argparse.ArgumentParser(description='1) Plot directions of your merged .h5 file;\n \
                                2) plot existing screen files w/o wcs frames; \n 3) find direction and its corresponding .h5 file name. \n \
                                Example 1) python3 Plot_dirandscreen.py --H5file h5file --size size --boxwidth boxwidth --FITSscreen fitsscreen \n \
                                Example 2) python3 Plot_dirandscreen.py --H5file h5file --size size --boxwidth boxwidth --FITSscreen fitsscreen --plot True \n \
                                Example 3) python3 Plot_dirandscreen.py --H5file h5file --H5file_ind merged_selfcalcyle \n')
#imaging options
parser.add_argument('--size', help='Size of the screen in pixels, default=64', default=64, type=int)
parser.add_argument('--boxwidth', help='Default=2.5, meaning the screen is covering 2.5x1.5 deg on each side. 1.5 is taken from Reinout\'s script, or \'ext_scale\' in this script', default=2.5, type=float)
parser.add_argument('--H5file', help='merged H5 file', type=str, required=True)
parser.add_argument('--H5file_ind', help='The first common series of letters of these individual H5 files, for example, merged_selfcalcyle', type=str)
parser.add_argument('--FITSscreen', help='Input FITS screen file; the RA and DEC would be taken from the screen .fits file', type=str)
parser.add_argument('--plot', help='Make screen plots for the input screen with RA and DEC', default=False)
#parser.add_argument('--ms', help='Measurement set to be imaged with IDG, phasecenter and antenna table are taken from the ms', type=str)
parser.add_argument('--plot_pix', help='Make screen plots for the input screen in pixels', default=False)
parser.add_argument('--RA', help='If you don\'t feed a screen using --FITSscreen, please type the RA of the field centre in decimal degree', default=0.0, type=float)
parser.add_argument('--DEC', help='If you don\'t feed a screen using --FITSscreen, please type the DEC of the field centre in decimal degree', default=0.0, type=float)


args        = vars(parser.parse_args())
h5file      =  args['H5file']
h5file_ind  = args['H5file_ind']
size        = args['size']   # default 64
boxwidth    = args['boxwidth']   #default 2.5 degree
#ms          = args['ms']
FITSscreen  = args['FITSscreen']
centre_RA   = args['RA']
centre_DEC  = args['DEC']
plot_wcs    = args['plot']
plot_pix    = args['plot_pix']
print (plot_wcs, plot_pix)

if FITSscreen is not None:
    hdul        = fits.open(FITSscreen)
    fits_data   = hdul[0].data
    fits_header = hdul[0].header
    centre_RA, centre_DEC = get_field_centre(fits_header)


print ('The current directory is ', PATH, 'in order to change the running directory, please change value of PATH in this script.')

antenna = get_antenna_names(h5file)


dir_array, dir_tag, dir_x, dir_y, dir_x_deg, dir_y_deg = get_direction(h5file)
dir_xy = np.column_stack((dir_tag, dir_x_deg, dir_y_deg))
dir_xy_pos = np.column_stack((dir_x, dir_y))
print ('There are ' + str(len(dir_array)) + ' directions.')
print ('These directions are:\n', dir_array)
print ('In degrees are:\n', dir_xy)

if plot_wcs == True:
    plot_screen_wcs(fits_data, fits_header, size, antenna)

if plot_pix == True:
    plot_screen(fits_data, fits_header, size, antenna)



print ('The direction and its corresponding .h5 file:')
dir_h5 = dir_h5_file(h5file_ind)
print (dir_h5)
print ('\n This list is saved in screen_plots/h5_names.txt.')

if centre_RA == 0 and centre_DEC == 0:
    raise Exception("No direction plot will be plotted, because the phase centre is not specified, you can either feed a .fits screen, or use --RA and --DEC to directly feed the phase centre in degrees")

plot_directions(h5file, centre_RA, centre_DEC, boxwidth)
#dir_x_pix, dir_y_pix, dir_tag_num = deg2pix(h5file, boxwidth, size)
#print (dir_x_pix, dir_y_pix, dir_tag_num)
