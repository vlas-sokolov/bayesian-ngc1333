#!/usr/bin/env python3
""" Make fits files with MLE parameters as cube slices """

from pyspecnest.chaincrunch import parcube as mle_parcube
from opencube import make_cube_shh

from config import name_id, npars, chain_dir, file_mle_formatter


def make_parcubes(peaks=(1, 2)):
    """ Generates parameter cubes with Maximum Likelihood Estimates """

    io_kwargs = dict(output_dir=chain_dir, name_id=name_id)

    cubes = make_cube_shh()

    mle_kwargs = dict(shape=cubes.cube.shape[1:],
                      npars=npars,
                      header=cubes.header)
    mle_kwargs.update(io_kwargs)

    for npeaks in peaks:
        mle_parcube(npeaks=npeaks, writeto=file_mle_formatter.format(npeaks),
                    mode='fast', stat='mle', **mle_kwargs)


if __name__ == '__main__':
    make_parcubes(peaks=(1, 2))
