#!/usr/bin/env python3
from pyspecnest.chaincrunch import cube_K, cube_Z
import opencube
import os

name_id = 'ngc1333-gas'
home_dir = os.path.expanduser('~')
proj_dir = os.path.join(home_dir, 'Projects', 'ngc1333-gas/')
chain_dir = os.path.join(proj_dir, 'nested-sampling' + os.path.sep)

lnKZ_xy_kwargs = dict(silent = True, output_dir = chain_dir,
                      name_id = name_id)

spc = opencube.make_cube()

# NOTE: why normalize is set to False? It's because it is more
# computationally expensive to include the normalization constant
# into likelihood calls, so on a normal sampling run we would get
# evidences from likelihoods that aren't normed.
cube_KZ_kwargs = dict(shape = spc.cube.shape[1:],
                      data = spc.cube,
                      peaks = [0, 1, 2],
                      rms = spc.errorcube,
                      header = spc.header,
                      normalize = False)

cube_KZ_kwargs.update(lnKZ_xy_kwargs)

hdu = cube_K(writeto=os.path.join(chain_dir, 'ngc1333-Ks.fits'), **cube_KZ_kwargs)
hdu = cube_Z(writeto=os.path.join(chain_dir, 'ngc1333-Zs.fits'), **cube_KZ_kwargs)
