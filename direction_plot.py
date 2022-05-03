# python3 

from astropy.io import fits
import tables
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
import os
import re
from astropy import units as u
from astropy.coordinates import Angle
from matplotlib.patches import Rectangle
from shapely.geometry.polygon import LinearRing, Polygon

#PATH = str(Path().absolute()).split("\/")[0] # Directory of current working directory
PATH = '/net/rijn2/data2/Haoyang/ALICE/myPybdsf/'
dir_path = PATH + '/plots/'
file_name = '/net/rijn2/data2/Haoyang/ALICE/myPybdsf/image_en1_field_1asec_facetallblocks_applyfacetbeam-MFS-image-pb.fits' 

if os.path.exists(dir_path)==False:
    os.mkdir(dir_path)

def get_image(image_name):
    hdul = fits.open(str(image_name))
    fits_data   = hdul[0].data
    image_data = fits_data[0,0,:,:]
    rms, rms0 = cal_RMS(image_data)
    return image_data, rms0
    
def deg2pix(direction_array_filename, wcs):
    '''
    convert degrees to pixels in the screen image plane using wcs function
    dir_x_pix: float list - dirction RAs in pixels
    dir_y_pix: float list - dirction DECs in pixels
    dir_tag_num: string list - direction tags only in number
    '''
    dir_x_pix = []
    dir_y_pix = []
    dir_array = np.load(str(direction_array_filename)) 
    for i in range(len(dir_array)):
        X= wcs.wcs_world2pix(dir_array[i][0], dir_array[i][1], 0)[0].tolist()
        Y= wcs.wcs_world2pix(dir_array[i][0], dir_array[i][1], 0)[1].tolist()
        dir_x_pix+=[X]
        dir_y_pix+=[Y]
    return dir_x_pix, dir_y_pix

def plot_polygon(direction_array_filename, facet_reg_file, centre_RA, centre_DEC, size, output_name, wcs):
    '''
    convert degrees to pixels in the screen image plane using wcs function
    dir_x_pix: float list - dirction RAs in pixels
    dir_y_pix: float list - dirction DECs in pixels
    dir_tag_num: string list - direction tags only in number
    '''
    fig = plt.figure(figsize=(5, 5))
    ax = plt.subplot(1,1,1, projection=wcs)
    dir_x_pix, dir_y_pix = deg2pix(direction_array_filename, wcs)
    centre_RA_pix, centre_DEC_pix = wcs.wcs_world2pix(centre_RA, centre_DEC, 0)
    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
    ax.grid(alpha=0.6)
    ax.scatter(centre_RA_pix, centre_DEC_pix, s = 10, edgecolor = 'k')
    plt.annotate('centre', (centre_RA_pix, centre_DEC_pix), ha='center', va='top', fontsize=14)
    ax.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 10, marker = 's', edgecolor = 'r')
    ax.set_xlim([-3500, 25000])
    ax.set_ylim([-3500, 25000])
    rect = Rectangle((0,0),size,size,linewidth=1,edgecolor='b',facecolor='none', alpha = 0.5)
    ax.add_patch(rect)
    ax.set_aspect(1)
    with open(facet_reg_file) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('polygon('):
            line_temp = re.findall(r'\d+\.\d+', line)
            line_num = [float(i) for i in line_temp]
            xy = [(x, y) for x, y in zip(line_num[::2], line_num[1::2])]
            poly_points = wcs.wcs_world2pix(xy,0)
            poly_coord = list(tuple(a.tolist()) for a in poly_points)
            poly = Polygon(poly_coord)
            x,y = poly.exterior.xy
            ax.plot(x, y, color='#6699cc', alpha=0.7, linewidth=1, solid_capstyle='round', zorder=2)
            #ax.set_title('Polygon')
    plt.savefig(dir_path + str(output_name))
    plt.show()

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

def plot_directions_wcs(centre_RA, centre_DEC, direction_array_filename, size, output_name):
    '''
    plot the directions within the wcs framework
    fits_header: the header of the .fits file, to take the phase centre values
    h5file: str, input merged h5parm name
    direction_array: numpy array, output of the get_direction function
    '''
    dir_x_pix, dir_y_pix = deg2pix(direction_array_filename, wcs)
    centre_RA_pix, centre_DEC_pix = wcs.wcs_world2pix(centre_RA, centre_DEC, 0)
    fig = plt.figure(figsize=(5, 5))
    ax = plt.subplot(1,1,1, projection=wcs)
    #ax.set_title('Direction Plot (wcs), boxwidth='+ str(boxwidth))
    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
    ax.grid(alpha=0.6)
    ax.scatter(centre_RA_pix, centre_DEC_pix, s = 10, edgecolor = 'k')
    #print (centre_RA_pix, centre_DEC_pix)
    #print (dir_x_pix, dir_y_pix)
    #print (max(dir_x_pix), min(dir_x_pix), max(dir_y_pix), min(dir_y_pix))
    print (size)
    plt.annotate('centre', (centre_RA_pix, centre_DEC_pix), ha='center', va='top', fontsize=14)
    ax.scatter(dir_x_pix, dir_y_pix, facecolors='none', s = 10, marker = 's', edgecolor = 'r')
    #for idx, txt in enumerate(dir_tag_num):
    #    ax.text(dir_x_pix[idx], dir_y_pix[idx], txt, style ='italic', color ="r")
    ax.set_xlim([-3500, 25000])
    ax.set_ylim([-3500, 25000])
    rect = Rectangle((0,0),size,size,linewidth=1,edgecolor='b',facecolor='none', alpha = 0.5)
    ax.add_patch(rect)
    ax.set_aspect(1)
    plt.savefig(dir_path + str(output_name))
    plt.show()
    
hdul = fits.open(file_name)
fits_data   = hdul[0].data
fits_header = hdul[0].header
wcs = WCS(fits_header, naxis = 2)
image_data = fits_data[0,0,:,:]
size = len(image_data)
centre_RA, centre_DEC = get_field_centre(fits_header)
print ('The phase centre in degrees are: ', str(centre_RA) + ', ' + str(centre_DEC))

facet_reg_file = '/net/rijn2/data2/Haoyang/ALICE/myPybdsf/facetsscreen.reg'
plot_directions_wcs(centre_RA, centre_DEC, PATH+'region_full.npy', size, 'full_dir.png')
#plot_directions_wcs(centre_RA, centre_DEC, PATH+'region_gy.npy', size, 'gy_dir.png')
plot_polygon(PATH+'region_gy.npy', facet_reg_file, centre_RA, centre_DEC, size, 'polygon.png', wcs)


