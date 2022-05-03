from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle
from regions import PixCoord
from regions import RectangleSkyRegion, RectanglePixelRegion
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import pyregion
from optparse import OptionParser
import argparse
from astropy.wcs import WCS
import numpy as np


def cut_recregion_pixel(x,y,width,height,angle,fits_data,file_name, figure_name):
  center = PixCoord(x=x, y=y)
  reg = RectanglePixelRegion(center=center, width=width,
			     height=height, angle=angle)
  mask = reg.to_mask()
  data = mask.cutout(fits_data)
  weighted_data = mask.multiply(fits_data)
  
  plt.subplots(1,3,figsize=(18,18)) 
  plt.subplot(1, 3, 1)
  plt.title("Mask")
  plt.imshow(mask.data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.subplot(1, 3, 2)
  plt.title("Data cutout")
  plt.imshow(data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.subplot(1, 3, 3)
  plt.title("Data cutout multiplied by mask")
  plt.imshow(weighted_data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.savefig(str(file_name) + '_demo.png', bbox_inches='tight')
  
  plt.figure()
  plt.title(str(figure_name))
  plt.imshow(data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)
  plt.xticks(color='w')
  plt.yticks(color='w')
  plt.savefig(str(file_name), bbox_inches='tight')
  
  return mask
  

def cut_recregion_reg(reg_file, fits_data, file_name, figure_name):
  regions = pyregion.open(str(reg_file))
  x = regions[0].coord_list[0]
  y = regions[0].coord_list[1]
  width = regions[0].coord_list[2]
  height = regions[0].coord_list[3]
  angle =  Angle(regions[0].coord_list[4], 'deg')
  center = PixCoord(x=x, y=y)
  reg = RectanglePixelRegion(center=center, width=width,
			     height=height, angle=angle)
  mask = reg.to_mask()
  data = mask.cutout(fits_data)
  weighted_data = mask.multiply(fits_data)
  
  '''
  fig, ax = plt.subplots(1, 1)
  patch = reg.as_artist(facecolor='none', edgecolor='red', lw=2)
  ax.add_patch(patch)

  plt.xlim(12300, 12800)
  plt.ylim(10000, 10500)
  ax.set_aspect('equal')
  plt.show()
  '''
  plt.subplots(1,3,figsize=(18,18)) 
  plt.subplot(1, 3, 1)
  plt.title("Mask")
  #plt.colorbar()
  plt.imshow(mask.data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.subplot(1, 3, 2)
  plt.title("Data cutout")
  plt.imshow(data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.subplot(1, 3, 3)
  plt.title("Data cutout multiplied by mask")
  plt.imshow(weighted_data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.savefig(str(file_name) + '_demo.png', bbox_inches='tight')
  
  plt.figure()
  plt.title(str(figure_name))
  plt.imshow(data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)
  plt.xticks(color='w')
  plt.yticks(color='w')
  plt.savefig(str(file_name), bbox_inches='tight')
  
  '''
  fig, ax = plt.subplots(1, 1, figsize=(12,12))
  ax.imshow(fits_data[int(x)-50:int(x)+50,int(y)-50:int(y)+50], cmap=plt.cm.viridis,
          interpolation='nearest', origin='lower')
  ax.add_artist(mask.bbox.as_artist(facecolor='yellow', edgecolor='white'))# outside of the image
  ax.add_artist(reg.as_artist(facecolor='red', edgecolor='orange'))
  plt.savefig(str(file_name) + '_background.png', bbox_inches='tight')
  '''
  return x,y,width,height,angle,mask


def cut_recregion_reg_wcs(reg_file, fits_data, fits_header, file_name, figure_name):
  regions = pyregion.open(str(reg_file))
  x = regions[0].coord_list[0]
  y = regions[0].coord_list[1]
  width = regions[0].coord_list[2]
  height = regions[0].coord_list[3]
  angle =  Angle(regions[0].coord_list[4], 'deg')
  center = PixCoord(x=x, y=y)
  reg = RectanglePixelRegion(center=center, width=width,
			     height=height, angle=angle)
  mask = reg.to_mask()
  data = mask.cutout(fits_data)
  np.savetxt(str(file_name)+'.csv', data, delimiter=",")
  weighted_data = mask.multiply(fits_data)
  wcs = WCS(fits_header, naxis = 2)
  
    
  plt.subplots(1,3,figsize=(18,18)) 
  plt.subplot(1, 3, 1)
  plt.title("Mask")
  #plt.colorbar()
  plt.imshow(mask.data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.subplot(1, 3, 2)
  plt.title("Data cutout")
  plt.imshow(data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.subplot(1, 3, 3)
  plt.title("Data cutout multiplied by mask")
  plt.imshow(weighted_data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)

  plt.savefig(str(file_name) + '_demo.png', bbox_inches='tight')
  
  fig = plt.figure()
  ax = plt.subplot(projection=wcs)
  ax.set_title(str(figure_name))
  ax.coords[1].set_ticklabel(exclude_overlapping=True)
  im = ax.imshow(data, cmap=plt.cm.viridis,
	     interpolation='nearest', origin='lower',
	     extent=mask.bbox.extent)
  ax.set_xlabel('Right Ascension (J2000)')
  ax.set_ylabel('Declination (J2000)')
  fig.colorbar(im, label=r'Flux ($\mu$Jy)')
  plt.savefig(str(file_name), bbox_inches='tight')
  
  
  return x,y,width,height,angle,mask

"""
hdul = fits.open('lockman_radio_image.fits')
hdu = hdul[0]

x, y = 12622, 10393
width, height = 15, 15
angle = Angle(0, 'deg')

mask = cut_recregion_pixel(x,y,width,height,angle,hdu.data[0,0,:,:], 'test1', 'LOFAR')  
mask = cut_recregion_reg('test2.reg',hdu.data[0,0,:,:], 'test2', 'LOFAR')
hdul.close()
"""

if __name__ == '__main__':
  string = '''Example: python Fits_cutout.py -o 'test2' --regfile 'test2.reg' --datafile 'lockman_radio_image.fits' --scopename "LOFAR" '''
  parser = OptionParser(usage=string)
  parser.add_option("-o", "--outputfile", type="str", dest="output",
		    help="Output png file name", metavar="FILE")
  parser.add_option("--regfile", type="str", dest="reg_filename",
		    help="Region (.reg) file name, should be saved in the format of physical not wcs", metavar="FILE")
  parser.add_option("--datafile", type="str", dest="data_filename",
		    help=".fits image file name", metavar="FILE")
  parser.add_option("--scopename", type="str", dest="scope_name",
		    help="Title of the cut-out image, can be the telescope/survey/source name", metavar="string", default = 'LOFAR')
  (options, args) = parser.parse_args()
  hdul = fits.open(options.data_filename)
  if hdul[0].data.all == None:
    hdu = hdul[1]
  else: 
    hdu = hdul[0]
  if len(hdu.data.shape) == 4:
    fits_input =  hdu.data[0,0,:,:]
  elif len(hdu.data.shape) == 2:
    fits_input =  hdu.data
  fits_header = hdu.header
  #x,y,width,height,angle,mask = cut_recregion_reg(options.reg_filename,fits_input, options.output, options.scope_name)
  x,y,width,height,angle,mask = cut_recregion_reg_wcs(options.reg_filename, fits_input, fits_header, options.output, options.scope_name)
 ## print (x,y,width,height,angle, mask.data[0], mask.data[0].shape)
  

 




