from astropy.io import fits
import tables
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
import os
#from pathlib import Path
from astropy.stats import sigma_clip
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle
from regions import PixCoord
from regions import RectangleSkyRegion, RectanglePixelRegion
from astropy.utils.data import get_pkg_data_filename
from matplotlib.patches import Rectangle
import pyregion
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


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

def crop_center(img,cropx,cropy):
    y,x = img.shape
    startx = x//2-(cropx//2)
    starty = y//2-(cropy//2)    
    return img[starty:starty+cropy,startx:startx+cropx]

def cal_RMS(image_data):
    #row_length, column_length = image_data.shape 
    #central_image = crop_center(image_data,row_length//10,column_length//10)
    rms = np.sqrt(np.mean(np.square(image_data)))
    #filtered_data = sigma_clip(image_data, sigma=5, maxiters=10)
    filtered_data = sigma_clip(image_data, sigma=5)
    rms0 = np.sqrt(np.mean(np.square(filtered_data)))
    print (rms, rms0)
    return rms, rms0


def plot_fits(fits_data, fits_header):
    '''
    plot the .fits image with wcs axises
    fits_data: numpy array taken from .fits file
    fits_header: the header of the .fits file, to take the phase centre values
    '''
    fig = plt.figure(figsize=(12, 9))
    image_data = fits_data[0,0,:,:]
    rms, rms0 = cal_RMS(image_data)
    wcs = WCS(fits_header, naxis = 2)
    #position = (size/2, size/2) # image centre
    ax = plt.subplot(1,1,1,projection=wcs)
    #ax.set_title('title')
    im = ax.imshow(1000*image_data, cmap=plt.cm.cubehelix,
             interpolation='nearest', origin='lower')
    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
    ax.grid(alpha=0.6)
    clb = plt.colorbar(im)
    clb.set_label(r'Flux density [mJy beam$^{-1}$]')
    im.set_clim(vmin=-1.5*rms0*1000, vmax=6*rms0*1000)
    #fig.colorbar(im)
    plt.savefig(dir_path + 'full_image.png', bbox_inches='tight', dpi = 300)
    plt.show()
    plt.close()


def plot_fits_centre(reg_file, fits_data, fits_header, output_name):
    regions = pyregion.open(str(reg_file))
    x = regions[0].coord_list[0]
    y = regions[0].coord_list[1]
    width = regions[0].coord_list[2]
    height = regions[0].coord_list[3]
    angle =  Angle(regions[0].coord_list[4], 'deg')
    center = PixCoord(x=x, y=y)
    reg = RectanglePixelRegion(center=center, width=width, height=height, angle=angle)
    mask = reg.to_mask()
    image_data = fits_data[0,0,:,:]
    data = mask.cutout(image_data)
    rms, rms0 = cal_RMS(data)
    wcs = WCS(fits_header, naxis = 2)
    fig = plt.figure(figsize=(12, 9))
    ax = plt.subplot(1,1,1,projection=wcs)
    #ax.set_title('title')
    im = ax.imshow(1000*data, cmap=plt.cm.cubehelix,
             interpolation='nearest', origin='lower', extent=mask.bbox.extent)
    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
    ax.grid(alpha=0.6)
    clb = plt.colorbar(im)
    clb.set_label(r'Flux density [mJy beam$^{-1}$]')
    im.set_clim(vmin=-1.5*rms0*1000, vmax=6*rms0*1000)
    #fig.colorbar(im)
    plt.savefig(dir_path + str(output_name), bbox_inches='tight', dpi = 300)
    plt.show()
    plt.close()


def plot_fits_reg(reg_file_list, im_name, output_name):
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(5, 8))
    image_data, rms = get_image(im_name)
    for ax, reg_file in zip(axes.flat, reg_file_list):
        ax.set_axis_off()
        regions = pyregion.open(str(reg_file))
        x = regions[0].coord_list[0]
        y = regions[0].coord_list[1]
        width = regions[0].coord_list[2]
        height = regions[0].coord_list[3]
        angle =  Angle(regions[0].coord_list[4], 'deg')
        center = PixCoord(x=x, y=y)
        reg = RectanglePixelRegion(center=center, width=width, height=height, angle=angle)
        mask = reg.to_mask()
        data = mask.cutout(image_data)
        rms, rms0 = cal_RMS(data)
        wcs = WCS(fits_header, naxis = 2)
        im = ax.imshow(1000*data, cmap='cubehelix',
                   vmin=-1.5*rms0*1000, vmax=12*rms0*1000)
        scalebar = AnchoredSizeBar(ax.transData, len(data)/10, '', 'lower left', 
                           pad=0.1,
                           color='white',
                           frameon=False,
                           size_vertical=1)
        ax.add_artist(scalebar)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), cax = fig.add_axes([0.78, 0.02, 0.03, 0.94]))
    cbar.set_label(r'Flux density [mJy beam$^{-1}$]')
    plt.subplots_adjust(left=0.40,
                    bottom=0.02, 
                    right=0.78, 
                    top=0.96, 
                    wspace=0, 
                    hspace=0.02)
    #plt.subplot_tool()
    #plt.title('')
    plt.show()
    plt.savefig(dir_path + str(output_name), dpi = 300)
    plt.close()
    
    
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

#plot_fits(fits_data, fits_header)
#reg_file = '/net/rijn2/data2/Haoyang/ALICE/myPybdsf/1arcsec/zoom.reg'
#plot_fits_centre(reg_file, fits_data, fits_header, output_name)
#plot_directions_wcs(centre_RA, centre_DEC, PATH+'region_full.npy', size, 'full_dir.png')
#plot_directions_wcs(centre_RA, centre_DEC, PATH+'region_gy.npy', size, 'gy_dir.png')

im_name1 = '/net/rijn2/data2/Haoyang/ALICE/myPybdsf/image_en1_field_1asec_facetallblocks_applyfacetbeam-MFS-image-pb.fits' 
im_name2 = '/net/rijn2/data2/Haoyang/ALICE/myPybdsf/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
im_name3 = '/net/rijn2/data2/Haoyang/ALICE/myPybdsf/en1_radio_image.fits'

reg_file_list = ['random1.reg', 'random2.reg', 'random3.reg', 'random4.reg']
reg_file_list2 = ['random1_6arcsec.reg', 'random2_6arcsec.reg', 'random3_6arcsec.reg', 'random4_6arcsec.reg']
reg_file_list3 = ['random1_deep.reg', 'random2_deep.reg', 'random3_deep.reg', 'random4_deep.reg']

plot_fits_reg(reg_file_list, im_name1, 'Comparison_1arcsec.png')
plot_fits_reg(reg_file_list2, im_name2, 'Comparison_6arcsec.png')
plot_fits_reg(reg_file_list3, im_name3, 'Comparison_deep.png')
