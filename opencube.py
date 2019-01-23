""" All spectral cube I/O helper functions reside here """

import os
import warnings
import pickle
import numpy as np
from astropy.io import fits
from astropy import log
import pyspeckit
from pyspeckit import cubes
from pyspeckit.spectrum.classes import units
from config import (file_nh311_dr1, file_nh322_dr1,
                    file_rms_nh311_dr1, file_rms_nh322_dr1, cube_save_kwargs)

def make_cube(files=(file_nh311_dr1, file_nh322_dr1),
              rms_files=(file_rms_nh311_dr1, file_rms_nh322_dr1)):
    """
    Opens the cube and calculates all the pre-fitting attributes of interest.
    """
    # make sure we're working on arrays (why?)
    files = np.atleast_1d([f for f in files])
    rms_files = np.atleast_1d([f for f in rms_files])

    if files.size > 1:
        spc_dict = {f: pyspeckit.Cube(f) for f in files}
        rmsmaps = {f: fits.getdata(ef) for f, ef in zip(files, rms_files)}
        for f in files:
            spc_dict[f].errorcube = np.repeat([rmsmaps[f]],
                                        spc_dict[f].xarr.size, axis=0)
        # now the errorcubes should merge automatically
        spc = pyspeckit.CubeStack([spc_dict[f] for f in files])
        spc.xarr.refX = spc.cubelist[0].xarr.refX
        spc.xarr.refX_unit = spc.cubelist[0].xarr.refX_unit
    else:
        spc = pyspeckit.Cube(files[0])
        rms = fits.getdata(rms_files[0])
        # easier to handle everything get_spectrum-related
        spc.errorcube = np.repeat([rms], spc.xarr.size, axis=0)

    # I don't see a reason why errorcube should be a masked array
    if type(spc.errorcube) == np.ma.MaskedArray:
        spc.errorcube = np.array(spc.errorcube)

    spc.xarr.velocity_convention = 'radio'
    spc.xarr.convert_to_unit('km/s')

    snr = (spc.cube / spc.errorcube).max(axis=0)

    # TODO: fix multinest-pipeline.py and run_multicube.py
    #spc.errmap = rms
    spc.snrmap = snr

    return spc


def make_cube_shh(**kwargs):
    """ Shush! Opens the cube without triggering a wall of warnings. """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        old_log = log.level
        log.setLevel('ERROR')
        spc = make_cube(**kwargs)
        log.setLevel(old_log)

    return spc


def save_xarr(xarr, target_dir=cube_save_kwargs['target_dir'], 
              target_xarr=cube_save_kwargs['target_xarr'],
              target_xarrkwargs=cube_save_kwargs['target_xarrkwargs'],
              saved_xarrkwargs=[
                  'unit', 'refX', 'refX_unit', 'center_frequency',
                  'center_frequency_unit', 'velocity_convention'],
              **kwargs):
    """
    Saves essential attributes of a SpectroscopicAxis instance to disk.
    The array is saves via `np.save`, and essential kwargs are pickled.

    Why do we need this?

    Opening a spectral cube every time an inference is made on a spectrum adds
    a lot of overheads. We have to write it in such a way that single spectra
    can be read in a lazy manner.

    Can't find a way to save `SpectroscopicAxes` in pyspeckit that doesn't
    break the irregularities from CubeStack's, and cubes don't have writers...
    """
    # not sure how to preserve SpectroscopicAxes instances
    if type(xarr) == units.SpectroscopicAxes:
        log.warning("Attempting to save a SpectroscopicAxes instance,"
                    " but the only way to read it will be as"
                    " a SpectroscopicAxis.")

    # save the actual xarr with all its irregularities, if any
    np.save(os.path.join(target_dir, target_xarr), xarr.value)

    # pickle all the informative xarr attributes
    xarrkwargs = {k:xarr.__getattribute__(k) for k in saved_xarrkwargs}
    with open(os.path.join(target_dir, target_xarrkwargs), 'wb') as f:
        pickle.dump(xarrkwargs, f)


def load_xarr(target_dir=cube_save_kwargs['target_dir'], 
              target_xarr=cube_save_kwargs['target_xarr'],
              target_xarrkwargs=cube_save_kwargs['target_xarrkwargs'],
              **kwargs):
    """
    Restores a Spectroscopic instance saved by `save_xarr` function
    """
    xarr_value = np.load(os.path.join(target_dir, target_xarr))

    with open(os.path.join(target_dir, target_xarrkwargs), 'rb') as f:
        xarrkwargs = pickle.load(f)

    xarr = units.SpectroscopicAxis(xarr_value, **xarrkwargs)

    return xarr


def save_datacube(spc, 
        target_dir=cube_save_kwargs['target_dir'],
        target_cubefile=cube_save_kwargs['target_cubefile'],
        target_errfile=cube_save_kwargs['target_errfile'], 
        target_header=cube_save_kwargs['target_header'],
        **kwargs):
    """
    In addition to SpectroscopicAxis, spectral cube data has to be accessed
    in a rapid manner to allow for efficient parallelization.

    Why aren't we storing a fits file instead again? Because dealing with
    irregular axes in it (e.g., ammonia (1,1)-(4,4)) lines together) is too
    problematic. In principle, multiple-HDU approach might be better, but
    it would be a better solution for pyspeckit CubeStack writers.
    """
    # save spc.xarr
    save_xarr(spc.xarr, target_dir=target_dir, **kwargs)

    # save the spectral cube data to a file
    # (the order of cubelst items is the same as in spc.xarr)
    np.save(os.path.join(target_dir, target_cubefile), spc.cube)
    # save the error cube data to a file
    # TODO: should I implicitly require this to be present?
    try:
        np.save(os.path.join(target_dir, target_errfile), spc.errorcube)
    except NotImplementedError: # really, numpy, really?
        assert type(spc.errorcube.data) == np.ndarray # errorcube is masked
        np.save(os.path.join(target_dir, target_errfile), spc.errorcube.data)

    header = spc.header.copy()
    header['CTARG'] = spc._first_cel_axis_num
    header['SYSTEM'] = spc.system
    with open(os.path.join(target_dir, target_header), 'wb') as f:
        pickle.dump(header, f)


def _header_cube_to_spectrum(h, x, y):
    """ Taken from SpectralCube.get_spectrum for consistency """

    ct = 'CTYPE{0}'.format(h['CTARG'])
    header = cubes.speccen_header(fits.Header(cards=[(k, v) for k, v in
                                              h.items() if k != 'HISTORY']),
                                  lon=x, lat=y, system=h['SYSTEM'],
                                  proj=(h[ct][-3:] if ct in h else 'CAR'))
    return header


def get_spectrum(x, y, 
                 target_dir=cube_save_kwargs['target_dir'],
                 target_xarr=cube_save_kwargs['target_xarr'],
                 target_xarrkwargs=cube_save_kwargs['target_xarrkwargs'],
                 target_cubefile=cube_save_kwargs['target_cubefile'],
                 target_errfile=cube_save_kwargs['target_errfile'], 
                 target_header=cube_save_kwargs['target_header'],
                 mmap_mode='r', **kwargs):
    """
    Fast initialization of (X, Y) spectra from a spectral cube.

    The header-making code block was shamelessly taken from pyspeckit,
    for consistency reasons.
    """
    xarr = load_xarr(target_dir, target_xarr, target_xarrkwargs)

    data = np.load(os.path.join(target_dir, target_cubefile),
                   mmap_mode=mmap_mode)[:, y, x]
    error = np.load(os.path.join(target_dir, target_errfile),
                    mmap_mode=mmap_mode)[:, y, x]

    with open(os.path.join(target_dir, target_header), 'rb') as f:
        h = pickle.load(f)
    header = _header_cube_to_spectrum(h, x, y)

    sp = pyspeckit.Spectrum(xarr=xarr, data=data, error=error, header=header)

    return sp


def clean_saved(target_dir=cube_save_kwargs['target_dir'],
                target_xarr=cube_save_kwargs['target_xarr'],
                target_xarrkwargs=cube_save_kwargs['target_xarrkwargs'],
                target_cubefile=cube_save_kwargs['target_cubefile'],
                target_errfile=cube_save_kwargs['target_errfile'], 
                target_header=cube_save_kwargs['target_header'],
                **kwargs):
    """
    Removes the xarr files used by `save_xarr` and `load_xarr` functions
    """
    os.remove(os.path.join(target_dir, target_xarr))
    os.remove(os.path.join(target_dir, target_xarrkwargs))
    os.remove(os.path.join(target_dir, target_cubefile))
    os.remove(os.path.join(target_dir, target_errfile))
    os.remove(os.path.join(target_dir, target_header))


def update_model(sp, fit_type='gaussian'):
    """
    Tie a model to a Cube/CubeStack. Should work for all the standard
    fitters; others can be added with Cube.add_fitter method.
    """
    try:
        allowed_fitters = sp.specfit.Registry.multifitters
        sp.specfit.fitter = allowed_fitters[fit_type]
    except KeyError:
        raise ValueError('Unsupported fit type: %s\n'
                         'Choose one from %s'
                         % (fit_type, allowed_fitters.keys()))
    log.info("Selected %s model" % fit_type)
    sp.specfit.fittype = fit_type
    sp.fittype = fit_type
