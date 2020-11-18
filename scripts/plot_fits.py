import numpy as np
import optparse
import os,sys
import re
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

o = optparse.OptionParser(usage='%prog [options] *.fits')
o.add_option('--size', dest='size', nargs=2, default=(8, 6), type="float", help='size of the output image (default: 8 6)')
o.add_option('--vmin', dest='vmin', default=None, type="float", help='minimum value of the colormap (default: data minimum)')
o.add_option('--vmax', dest='vmax', default=None, type="float", help='maximum value of the colormap (default: data maximum)')
o.add_option('--cmap', dest='cmap', default="jet", type="str", help='colormap (default: jet)')
o.add_option('--cbar_label', dest='cbar_label', default=None, type="str", help='label of the color bar (default: BUNIT from header)')
o.add_option('--gridc', dest='gridc', default="silver", type="str", help='coordinates grid color (default: silver)')
o.add_option('--title', dest='title', default=None, type="str", help='title of the image (default: fits name)')
g = optparse.OptionGroup(o, 'Saving Options')
g.add_option('--pdf', dest='save_pdf', default=False, action='store_true', help='save image as pdf')
g.add_option('--png', dest='save_png', default=False, action='store_true', help='save image as png')
g.add_option('--name', dest='outfile', default=None, type="str", help='output image name with no extension (default: fits name)')
o.add_option_group(g)
opts,args = o.parse_args(sys.argv[1:])

fitsname = args[0]
hdul = fits.open(fitsname)
hdr = hdul[0].header
data = hdul[0].data[0,0,:,:]
wcs = WCS(hdr,naxis=2)

fig = plt.figure(figsize=opts.size)
ax = fig.add_subplot(1,1,1, projection=wcs)
ra = ax.coords[0]
dec = ax.coords[1]

if opts.vmin == None: opts.vmin = np.min(data)
if opts.vmax == None: opts.vmax = np.max(data)
if opts.cbar_label == None: opts.cbar_label = hdr['BUNIT']

plt.imshow(data, origin='lower', cmap=opts.cmap,
	aspect='auto', interpolation='none', 
	vmin=opts.vmin, vmax=opts.vmax)
plt.colorbar().set_label(opts.cbar_label, rotation=90)
ax.coords.grid(color=opts.gridc, ls='solid', lw=0.50)
plt.xlabel('Right Ascension (J2000)')
plt.ylabel('Declination (J2000)')

if opts.outfile == None: opts.outfile = re.split('/|.fits', fitsname)[-2]

if opts.title == None: opts.title = opts.outfile
plt.title(opts.title)
if opts.save_pdf:
	plt.savefig(opts.outfile + '.pdf')
if opts.save_png:
	plt.savefig(opts.outfile + '.png', dpi=200)
plt.show()
plt.close()

