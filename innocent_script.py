"""
The main sampling script.

Takes pixel positions (x, y) and the number of components (npeaks)
from the command line (via sys.argv[]), sets the priors, the likelihood
function, preforms some sanity checks, and fires up MultiNest through
pymultinest.
"""

import os
import sys
import warnings
import mpi4py
import numpy as np
from astropy import log
from astropy.io import fits
# Those import warnings are annoying
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    import pymultinest
    from pyspeckit.spectrum.models.ammonia import cold_ammonia_model
    from pyspecnest.ammonia import get_nh3_model
    from pyspecnest.chaincrunch import pars_xy, lnZ_xy, get_zero_evidence

# All the I/O functions now reside here
import opencube
# Configuration for line modelling setup (gets passed to pyspeckit)
from config import line_names, npars
# Path settings
from config import name_id, proj_dir, file_Zs
from config import default_yx, default_npeaks
# Kwargs for digesing and saving spectral cube (making it on the
# fly every time we need a spectrum is too slow!)
from config import cube_save_kwargs
# Finally, MultiNest settings and priors
from config import n_live_points, sampling_efficiency, get_priors_xoff_wrapped

# Compatibility with Python 2.7
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def mpi_rank():
    """
    Returns the rank of the calling process.
    """
    comm = mpi4py.MPI.COMM_WORLD
    rank = comm.Get_rank()

    return rank

# pop culture references always deserve their own function
def i_am_root():
    """
    Checks if the running subprocess is of rank 0
    """
    try:
        return True if mpi_rank() == 0 else False
    except AttributeError:
        # not running MPI
        return True

try:
    npeaks = int(sys.argv[1])
except IndexError:
    npeaks = default_npeaks
    log.info("npeaks not specified, setting to {}".format(npeaks))

try:
    yx = int(sys.argv[2]), int(sys.argv[3])
except IndexError:
    yx = default_yx
    log.info("xy-pixel not specified, setting to {}".format(yx[::-1]))

try:
    plotting = bool(int(sys.argv[4]))
except IndexError:
    plotting = 0

if not plotting:
    plot_fit, plot_corner, show_fit, show_corner = False, False, False, False
else: # defaults a for non-batch run
    plot_fit = True
    plot_corner = True
    show_fit = True
    show_corner = True
    from chainconsumer import ChainConsumer  # Optional if no plotting is done
    import matplotlib.pyplot as plt
    plt.rc('text', usetex=True)

y, x = yx
sp = opencube.get_spectrum(x, y, **cube_save_kwargs)

fittype_fmt = 'cold_ammonia_x{}'
fitmodel = cold_ammonia_model
sp.specfit.Registry.add_fitter(fittype_fmt.format(npeaks), npars=npars,
                               function=fitmodel(line_names=line_names))
# is this still needed?
opencube.update_model(sp, fittype_fmt.format(npeaks))

# needed because of this:
# https://github.com/pyspeckit/pyspeckit/issues/179
sp.specfit.fitter.npeaks = npeaks
# npeaks > 1 seems to break because of fitter.parnames is not mirrored
if len(sp.specfit.fitter.parnames) == npars and npeaks > 1:
    sp.specfit.fitter.parnames *= npeaks

# TODO: make a wrapper function for pyspecnest instead!
priors = get_priors_xoff_wrapped(npeaks)

nh3_model = get_nh3_model(sp, line_names, sp.error,
                          priors=priors, npeaks=npeaks)

# Safeguard - check some common causes of failure before scheduling
# a job that would just throw tens of thousands of errors at us
no_valid_chans = not np.any(np.isfinite(nh3_model.ydata))
sanity_check = np.isfinite(
    nh3_model.log_likelihood([15, 5, 15, 0.2, 7, 0.5] * npeaks,
                             nh3_model.npars, nh3_model.dof))
if no_valid_chans or not sanity_check:
    # This should fail if, e.g., the errors are not finite
    log.error("no valid pixels at x={}; y={}. Aborting.".format(*yx[::-1]))
    sys.exit()

# The first process gets to make the directory structure!
output_dir = os.path.join(proj_dir, 'nested-sampling/')
fig_dir = os.path.join(output_dir, 'figs/')
suffix = 'x{1}y{0}'.format(*yx)
chains_dir = '{}chains/{}_{}'.format(output_dir, name_id, suffix)
if not os.path.exists(chains_dir):
    try: # hacks around a race condition
        os.makedirs(chains_dir)
    except OSError as e:
        if e.errno != 17:
            raise

chains_dir = '{}/{}-'.format(chains_dir, npeaks)

# Run MultiNest on the model+priors specified
pymultinest.run(nh3_model.xoff_symmetric_log_likelihood,
                nh3_model.prior_uniform, nh3_model.npars,
                outputfiles_basename=chains_dir,
                verbose=True, n_live_points=n_live_points,
                sampling_efficiency=sampling_efficiency)

# The remainder of the script is not essential for sampling, and can be safely
# moved out into a script of its own.
if i_am_root() and plot_fit:
    # parse the results as sensible output
    from pyspecnest.chaincrunch import analyzer_xy
    a = analyzer_xy(x, y, npeaks, output_dir=output_dir,
                    name_id=name_id, npars=npars)

    a_lnZ = a.get_stats()['global evidence']
    log.info('ln(Z) for model with {} line(s) = {:.1f}'.format(npeaks, a_lnZ))

    try:
        lnZ0 = fits.getdata(file_Zs)[0]
    except (FileNotFoundError, OSError) as e:
        cubes = opencube.make_cube_shh()
        lnZ0 = get_zero_evidence(data=cubes.cube, rms=cubes.errorcube,
                                 normalize=False)

    Zs = lnZ_xy(list(np.arange(npeaks+1)), x=x, y=y, output_dir=output_dir,
                name_id=name_id, silent=True, lnZ0=(lnZ0[y, x], 0))

    log.info('ln(Z{}/Z{}) = {:.2f}'.format(npeaks, npeaks-1,
                                           Zs[npeaks] - Zs[npeaks-1]))
    if npeaks > 1:
        log.info('ln(Z{}/Z{}) = {:.2f}'.format(npeaks, 0, Zs[npeaks] - Zs[0]))

if plot_fit and i_am_root():
    sp.plotter(errstyle='fill')

    mle_pars = pars_xy(x=x, y=y, npars=npars, npeaks=npeaks,
                       output_dir=output_dir, name_id=name_id)

    mle_parinfo = sp.specfit.fitter._make_parinfo(mle_pars, npeaks=npeaks)[0]

    try:
        sp.specfit.plot_fit(xarr=sp.xarr, pars=mle_parinfo,
                            show_components=True)
    except TypeError:
        # eh? does it want pars or parinfo?
        sp.specfit.plot_fit(xarr=sp.xarr, pars=mle_pars, show_components=True)

    # annotate the Bayes factors
    plt.annotate('ln(Z{}/Z{}) = {:.2f}'.format(npeaks, npeaks-1,
                 Zs[npeaks] - Zs[npeaks-1]), xy=(0.05, 0.90),
                 xycoords='axes fraction')
    if npeaks > 1:
        plt.annotate('ln(Z{}/Z{}) = {:.2f}'.format(npeaks, 0,
                     Zs[npeaks] - Zs[0]), xy=(0.05, 0.85),
                     xycoords='axes fraction')
    if show_fit:
        plt.show()

    fig_name = "{}-fit-{}-x{}".format(name_id, suffix, npeaks)
    plt.savefig(os.path.join(fig_dir, fig_name + ".pdf"))

if plot_corner and i_am_root():
    mle_multinest = pars_xy(x=x, y=y, npars=npars, npeaks=npeaks,
                            output_dir=output_dir, name_id=name_id)
    unfrozen_slice = nh3_model.get_nonfixed_slice(a.data.shape, axis=1)
    c = ChainConsumer()
    parameters = nh3_model.get_names(latex=True, no_fixed=True)
    c.add_chain(a.data[:, 2:][unfrozen_slice], parameters=parameters)
    c.configure(statistics="max", summary=True)
    fig = c.plotter.plot(figsize="column")
    fig.get_size_inches()
    fig.set_size_inches(9, 7)

    fig_name = "{}-corner-{}-x{}".format(name_id, suffix, npeaks)
    plt.savefig(fig_dir + fig_name + ".pdf")

    if show_corner:
        plt.show()
