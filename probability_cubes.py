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

def write_2dhist(fname, idx, nbins, vmin, vmax, ctype3='', cunit3=''):
    def histogramize_post(arr):
        return np.histogram(arr[:, idx+2], density=True, bins=nbins,
                            range=(vmin, vmax))

    post_map = np.full(shape=((nbins,)+sort_arr.shape), fill_value=np.nan)
    histres = cubes.parallel_map(histogramize_post, resarray,
                                 numcores=multicore)
    for ((x, y), h) in zip(xylist, histres):
        post_map[:, y, x] = h[0]

    normed_post_map = post_map / np.nansum(post_map, axis=0)

    head = fits.getheader(bfactor_fname)
    head['CTYPE3'] = ctype3
    head['CUNIT3'] = cunit3
    head['NAXIS3'] = nbins
    head['BUNIT'] = 'PROBABILITY MASS FUNCTION'
    head['CRPIX3'] = 1
    head['CRVAL3'] = vmin
    head['CDELT3'] = (vmax - vmin) / nbins
    for key in ['PLANE1', 'PLANE2', 'PLANE3', 'PLANE4']:
        head.pop(key)

    hdu = fits.PrimaryHDU(data=normed_post_map, header=head)
    hdu.writeto(fname, overwrite=True)

suff = {1: [''],
        2: ['_1', '_2']}
for i, s in enumerate(suff[npeaks]):
    fname_fmt = 'probability-cubes/{}_x{}' + s + '.fits'
    for fname, idx, nbins, vmin, vmax, ctype3, cunit3 in zip(
            [fname_fmt.format(key, npeaks)
                for key in ['tkin', 'ntot', 'sig', 'vlsr']],
            [0, 2, 3, 4], [500]*4,
            [3 , 12, 0.08, 5.1],
            [24, 15, 1.00, 9.9],
            ['KINETIC TEMPERATURE', 'COLUMN DENSITY',
                'VELOCITY DISPERSION', 'CENTROID VELOCITY'],
            ['K', 'cm-2', 'km s-1', 'km s-1']):
        write_2dhist(fname, idx+i*npars, nbins, vmin, vmax, ctype3, cunit3)
