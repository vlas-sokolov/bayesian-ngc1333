import os
import numpy as np
from astropy import log
from astropy.io import fits
from pyspeckit import cubes
from pyspecnest.chaincrunch import analyzer_xy
from pyspecnest.pool_multinest import get_xy_sorted

npars = 6
npeaks = 1
name_id = 'ngc1333-gas'
proj_dir = os.path.expanduser('~/Projects/ngc1333-gas/')
output_dir = os.path.join(proj_dir, 'nested-sampling/')
bfactor_fname = os.path.join(output_dir, 'ngc1333-Ks.fits')
multicore = 30

sort_arr = fits.getdata(bfactor_fname)[0]
xylist = get_xy_sorted(sort_arr)

def get_xy_data(args):
    x, y, a = args
    return x, y, a.get_data()

post_file = os.path.join(output_dir,
                '{}-x{}-posterior.npy'.format(name_id, npeaks))
try:
    resarray = np.load(post_file)
except (FileNotFoundError, OSError) as e:
    ankwargs = dict(npeaks=npeaks, output_dir=output_dir,
                    name_id=name_id, npars=npars)
    xy_analyzer_list = [(x, y, analyzer_xy(x, y, **ankwargs))
                        for (x, y) in xylist]

    # runs for ~5-10 minutes
    result = cubes.parallel_map(get_xy_data, xy_analyzer_list,
                                numcores=multicore)
    resarray = np.array([arr for (x, y, arr) in result], dtype=object)
    # takes a while to save as well
    np.save(post_file, resarray)

nbins = 500
vlsr_min, vlsr_max = 5.1, 9.9
def histogramize(arr):
    return np.histogram(arr[:, 6], density=True, bins=nbins,
                        range=(vlsr_min, vlsr_max))

# very manageble
vlsr_map = np.full(shape=((nbins,)+sort_arr.shape), fill_value=np.nan)
histres = cubes.parallel_map(histogramize, resarray, numcores=multicore)
for ((x, y), h) in zip(xylist, histres):
    vlsr_map[:, y, x] = h[0]

normed_vlsr_map = vlsr_map / np.nansum(vlsr_map, axis=0)

head = fits.getheader(bfactor_fname)
head['CTYPE3'] = 'VELOCITY'
head['CUNIT3'] = 'km/s'
head['NAXIS3'] = nbins
head['BUNIT'] = 'PROBABILITY MASS FUNCTION'
head['CRPIX3'] = 1
head['CRVAL3'] = vlsr_min
head['CDELT3'] = (vlsr_max - vlsr_min) / nbins
for key in ['PLANE1', 'PLANE2', 'PLANE3', 'PLANE4']:
    head.pop(key)

# save to file
hdu = fits.PrimaryHDU(data=normed_vlsr_map, header=head)
hdu.writeto('vlsr_test.fits', overwrite=True)


# now let's do the same for kinetic temperature
nbins = 500
tkin_min, tkin_max = 2.725, 25
def histogramize_tkin(arr):
    return np.histogram(arr[:, 2], density=True, bins=nbins,
                        range=(tkin_min, tkin_max))

tkin_map = np.full(shape=((nbins,)+sort_arr.shape), fill_value=np.nan)
histres = cubes.parallel_map(histogramize_tkin, resarray, numcores=multicore)
for ((x, y), h) in zip(xylist, histres):
    tkin_map[:, y, x] = h[0]

normed_tkin_map = tkin_map / np.nansum(tkin_map, axis=0)

head = fits.getheader(bfactor_fname)
head['CTYPE3'] = 'TEMPERATURE'
head['CUNIT3'] = 'K'
head['NAXIS3'] = nbins
head['BUNIT'] = 'PROBABILITY MASS FUNCTION'
head['CRPIX3'] = 1
head['CRVAL3'] = tkin_min
head['CDELT3'] = (tkin_max - tkin_min) / nbins
for key in ['PLANE1', 'PLANE2', 'PLANE3', 'PLANE4']:
    head.pop(key)

# testing
hdu = fits.PrimaryHDU(data=normed_tkin_map, header=head)
hdu.writeto('tkin_test.fits', overwrite=True)
