#!/usr/bin/env python3
from pyspecnest.chaincrunch import parcube as mle_parcube
import os
from opencube import make_cube_shh

name_id = 'ngc1333-gas'
home_dir = os.path.expanduser('~')
proj_dir = os.path.join(home_dir, 'Projects', 'ngc1333-gas/')
chain_dir = os.path.join(proj_dir, 'nested-sampling' + os.path.sep)
npars = 6

io_kwargs = dict(output_dir=chain_dir, name_id=name_id)

cubes = make_cube_shh()

mle_kwargs = dict(shape = cubes.cube.shape[1:],
                  npars = npars,
                  header = cubes.header)
mle_kwargs.update(io_kwargs)

writeto = os.path.join(chain_dir, '{}-mle'.format(name_id) + '-x{}.fits')

peaks = [1, 2]
for npeaks in peaks:
    mle_parcube(npeaks=npeaks, writeto=writeto.format(npeaks),
                mode='fast', stat='mle', **mle_kwargs)
