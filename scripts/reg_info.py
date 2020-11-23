import optparse
import os,sys
import regions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

o = optparse.OptionParser(usage='%prog [options] *.fits')
o.add_option('-r', '--reg', dest='reg', type="str", help='region file (allowed extensions: reg, crtf, fits)')
o.add_option('-p', '--plot', dest='plot', action="store_true", help='simple plot of the selected region')
opts,args = o.parse_args(sys.argv[1:])

for fitsname in args:
    hdul = fits.open(fitsname)
    hdr = hdul[0].header
    data = hdul[0].data[0,0,:,:]
    wcs = WCS(hdr,naxis=2)
    fits_units = hdr['BUNIT']
    hdul.close()

    reg_type = opts.reg.split('.')[-1]

    if reg_type == 'reg': read_reg = regions.read_ds9
    elif reg_type == 'crtf': read_reg = regions.read_crtf
    elif reg_type == 'fits': read_reg = regions.read_fits_region
    else: raise Exception('invalid region type (only reg, crtf and fits extensions are accepted)')

    all_regions = read_reg(opts.reg)
    reg_sum = []; reg_rms = []; reg_mean = []; reg_ra = []; reg_dec = []; reg_npoints = []
    idx = 0
    for region in all_regions:
        reg_pxl = region.to_pixel(wcs=wcs)
        mask = reg_pxl.to_mask()
        data_reg = mask.cutout(data)
        if opts.plot:
            plt.imshow(data_reg, cmap='viridis', origin='lower')
            plt.colorbar(label=fits_units)
            plt.title('Region ' + str(idx))
            plt.axis('off')
            plt.show()
            idx += 1
        reg_ra = np.append(reg_ra, str(SkyCoord(region.center).ra.to('hourangle')))
        reg_dec = np.append(reg_dec, str(SkyCoord(region.center).dec))  
        reg_npoints = np.append(reg_npoints, data_reg.size)
        reg_sum = np.append(reg_sum, np.sum(data_reg))
        reg_mean = np.append(reg_mean, np.mean(data_reg))
        reg_rms = np.append(reg_rms, np.std(data_reg))
 
    reg_table = {   'central RA': reg_ra, 
                    'central dec': reg_dec, 
                    'npoints': reg_npoints, 
                    'sum ['+fits_units+']' : reg_sum, 
                    'mean ['+fits_units+']' : reg_mean, 
                    'rms ['+fits_units+']' : reg_rms 
                    }

    reg_df = pd.DataFrame(reg_table)
    reg_df.index.name = 'Region'

    print('\n=== Image ' + fitsname + ' with region file ' + opts.reg + ' ===\n')
    print(reg_df)
    print('-'*reg_df.size)
    print('Total points:\t', np.sum(reg_npoints))
    print('Total sum:\t%.2e'%np.sum(reg_sum), fits_units)
    print('Average sum:\t%.2e'%np.mean(reg_sum), fits_units)
    print('Average mean:\t%.2e'%np.mean(reg_mean), fits_units)
    print('Average rms:\t%.2e'%np.mean(reg_rms), fits_units)
