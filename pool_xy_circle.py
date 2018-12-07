"""
Instead of MPI parallelisation for each pixel, it is [citation needed] faster
to simply distribute pixels to different processes via pooling. The speedup is
not about the actual sampling, but the overheads are only executed once...
"""
from __future__ import division
# the os.niceness will be inherited by child processes
# behave, kids!
import os
os.nice(19)
import sys
import multiprocessing
from collections import OrderedDict
import numpy as np
from pyspecnest import pool_multinest
from opencube import make_cube_shh
from config import sampler_script_file


def try_get_args(n, fallback, forcetype=str):
    try:
        # sys.argv[0] is some env executable path...
        arg = forcetype(sys.argv[n+1])
    except IndexError:
        arg = fallback

    return arg

def get_circle_mask(d, shape=None, x0=None, y0=None):
    if shape is None:
        shape = (d, d) # will make a square
    if y0 is None:
        y0 = d / 2.
    if x0 is None:
        x0 = d / 2.

    ymax, xmax = shape
    rad_yy, rad_xx = np.meshgrid(np.arange(xmax) - x0 + .5,
                                 np.arange(ymax) - y0 + .5)
    dist = (rad_yy**2 + rad_xx**2)**0.5
    dist[dist > (d / 2.)] = np.nan
    return np.isfinite(dist)

def get_circle_idf(x0, y0, d):
    """ Returns x- and y-indices of a circe in an array """
    idx, idy = np.where(get_circle_mask(d))

    return idx + x0 - d / 2., idy + y0 - d / 2.

def main():
    # NOTE: normal dict would mess up the order of the arguments
    default_args = OrderedDict([('npeaks', 2), ('method', "snr"), ('cut', 8),
                                ('n_cpu', 15), ('x0', 113), ('y0', 140), ('r', 60)])

    runtime_args = {}
    for i, (argname, argval) in enumerate(default_args.items()):
        runtime_args[argname] = try_get_args(i, argval, type(argval))

    method = runtime_args.pop('method')
    n_cpu = runtime_args.pop('n_cpu')
    npeaks = runtime_args['npeaks']
    cut = runtime_args['cut']
    x0, y0 = runtime_args['x0'], runtime_args['y0']
    r = runtime_args['r']

    spc = make_cube_shh() # comes with pregen snr attributes...
    mask = get_circle_mask(r*2, spc.snrmap.shape, x0, y0)
    masked_snr = spc.snrmap.copy()

    # I reap what I sow - a result of over-hacking things up
    assert cut >= 1

    masked_snr[~mask] = cut - 1
    order = pool_multinest.get_xy_sorted(masked_snr, cut=cut)
    tasklist_kwargs=dict(n_cpu=n_cpu, method='snr')

    tasklist_kwargs.update(runtime_args)

    tasks = pool_multinest.get_tasks(n_cpu, xy_order=order,
            npeaks=npeaks, script=sampler_script_file)

    pool = multiprocessing.Pool(processes=n_cpu)
    # NOTE: map vs imap? imap has better ordering... see more here:
    #       [https://stackoverflow.com/questions/26520781]
    # NOTE 2: imap won't work here...
    pool.map(pool_multinest.work, tasks)


if __name__ == '__main__':
    main()
