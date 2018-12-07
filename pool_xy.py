"""
Instead of MPI parallelisation for each pixel, it is [citation needed] faster
to simply distribute pixels to different processes via pooling. The speedup is
not about the actual sampling, but the overheads are only executed once...
"""
from __future__ import division
import sys
# the os.niceness will be inherited by child processes
# behave, kids!
import os
os.nice(19)
import multiprocessing
from collections import OrderedDict
import numpy as np
from pyspecnest import pool_multinest
from opencube import make_cube
from config import file_Ks, sampler_script_file

def try_get_args(n, fallback, forcetype=str):
    try:
        # sys.argv[0] is some env executable path...
        arg = forcetype(sys.argv[n+1])
    except IndexError:
        arg = fallback

    return arg

def main():
    # NOTE: normal dict would mess up the order of the arguments
    default_args = OrderedDict([('npeaks', 1), ('method', 'snr'), ('cut', 8),
                                ('n_cpu', 7)])

    runtime_args = {}
    for i, (argname, argval) in enumerate(default_args.items()):
        runtime_args[argname] = try_get_args(i, argval, type(argval))

    method = runtime_args.pop('method')
    n_cpu = runtime_args.pop('n_cpu')
    npeaks = runtime_args['npeaks']
    cut = runtime_args['cut']

    if method == 'snr':
        spc = make_cube() # comes with pregen snr attributes...
        sort_array = spc.snrmap
    elif method == 'Bfactor':
        from astropy.io import fits
        sort_array = fits.getdata(file_Ks)[npeaks-2]
    elif method == 'chisq':
        raise NotImplementedError

    order = pool_multinest.get_xy_sorted(sort_array,
                np.indices(sort_array.shape), cut=cut)

    tasks = pool_multinest.get_tasks(n_cpu, xy_order=order,
            npeaks=npeaks, script=sampler_script_file)

    pool = multiprocessing.Pool(processes=n_cpu)
    # NOTE: map vs imap? imap has better ordering... see more here:
    #       [https://stackoverflow.com/questions/26520781]
    # NOTE 2: imap won't work here...
    pool.map(pool_multinest.work, tasks)

if __name__ == '__main__':
    main()
