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

ankwargs = dict(npeaks=npeaks, output_dir=output_dir,
                name_id=name_id, npars=npars)
xy_analyzer_list = [(x, y, analyzer_xy(x, y, **ankwargs)) for (x, y) in xylist]

def get_xy_data(args):
    x, y, a = args
    return x, y, a.get_data()

post_file = os.path.join(output_dir,
                '{}-x{}-posterior.npy'.format(name_id, npeaks))
try:
    resarray = np.load(post_file)
except (FileNotFoundError, OSError) as e:
    # runs for ~5-10 minutes
    result = cubes.parallel_map(get_xy_data, xy_analyzer_list,
                                numcores=multicore)
    resarray = np.array([arr for (x, y, arr) in result], dtype=object)
    # takes a while to save as well
    np.save(post_file, resarray)

# brute-force way: make a big-ass array and put everything there
#a = np.full(shape=((int(np.nanmax(chain_lengths)),) + sort_arr.shape), fill_value=np.nan)
#for ((x, y), r) in zip(xylist, resarray):
#    a[:r[:, 4].size, y, x] = r[:, 4]

#bins = np.linspace(5.1, 9.9, 501)
def histogramize(arr):
    return np.histogram(arr[:, 6], density=True, bins=500, range=(5.1, 9.9))

# very manageble
vlsr_map = np.full(shape=((500,)+sort_arr.shape), fill_value=np.nan)
histres = cubes.parallel_map(histogramize, resarray, numcores=multicore)
for ((x, y), h) in zip(xylist, histres):
    vlsr_map[:, y, x] = h[0]

normed_vlsr_map = vlsr_map / np.nansum(vlsr_map, axis=0)

head = fits.getheader(bfactor_fname)
head['NAXIS3'] = 500
head['BUNIT'] = 'km/s'
head['CRPIX3'] = 1
head['CRVAL3'] = np.linspace(5.1, 9.9, 501)[0]
head['CDELT3'] = np.diff(np.linspace(5.1, 9.9, 501))[0]
for key in ['PLANE1', 'PLANE2', 'PLANE3', 'PLANE4']:
    head.pop(key)

# testing
hdu = fits.PrimaryHDU(data=normed_vlsr_map, header=head)
hdu.writeto('vlsr_test.fits', overwrite=True)
