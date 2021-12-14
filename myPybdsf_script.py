import bdsf

# Frits: pass the PB-corrected image to fitsname and the not-PB-corrected image to detection_image
fitsname = 'image_en1_field_1asec_screen-MFS-image-pb.fits'
detectimage = 'image_en1_field_1asec_screen-MFS-image.fits'

# Obtain the frequency from the fits header, CRVAL3
restfrq = 144627380.371094

img = bdsf.process_image(fitsname, detection_image=detectimage,
			thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', 
			ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15),
			group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, 
			atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,
			advanced_opts=True, blank_limit=None, frequency=restfrq)

#img.show_fit(ch0_islands = True, ch0_image = True)
# mask of islands (0 = outside island, 1 = inside island)
img.export_image(outfile = 'island_mask.fits', img_format = 'fits', img_type = 'island_mask')
# rms map image
img.export_image(outfile = 'rms.fits', img_format = 'fits', img_type = 'rms')
# mean map image
img.export_image(outfile = 'mean.fits', img_format = 'fits', img_type = 'mean')
# Gaussian model residual image
img.export_image(outfile = 'gaus_resid.fits', img_format = 'fits', img_type = 'gaus_resid')
# Gaussian model image
img.export_image(outfile = 'gaus_model.fits', img_format = 'fits', img_type = 'gaus_model')
# image used for source detection
img.export_image(outfile = 'ch0.fits', img_format = 'fits', img_type = 'ch0')

img.write_catalog(outfile = 'catalog.fits', format = 'fits', catalog_type =  'srl')
img.write_catalog(outfile = 'catalog.ascii', format = 'ascii', catalog_type =  'srl')
img.write_catalog(outfile = 'catalog.csv', format = 'csv', catalog_type =  'srl')

img.write_catalog(outfile = 'catalog_gau.fits', format = 'fits', catalog_type = 'gaul')
img.write_catalog(outfile = 'catalog_gau.ascii', format = 'ascii', catalog_type = 'gaul')
img.write_catalog(outfile = 'catalog_gau.csv', format = 'csv', catalog_type = 'gaul')
