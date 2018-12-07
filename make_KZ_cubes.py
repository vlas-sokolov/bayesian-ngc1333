#!/usr/bin/env python3
""" Make fits files with Bayesian evidences and Bayes factors """

from pyspecnest.chaincrunch import cube_K, cube_Z
import opencube

from config import name_id, chain_dir, file_Ks, file_Zs


def make_KZ_cubes(peaks=(0, 1, 2)):
    """ Generates evidence and Bayes factor maps """

    lnKZ_xy_kwargs = dict(silent=True, output_dir=chain_dir,
                          name_id=name_id)

    spc = opencube.make_cube()

    # NOTE: why normalize is set to False? It's because it is more
    # computationally expensive to include the normalization constant
    # into likelihood calls, so on a normal sampling run we would get
    # evidences from likelihoods that aren't normed.
    cube_KZ_kwargs = dict(shape=spc.cube.shape[1:],
                          data=spc.cube,
                          peaks=peaks,
                          rms=spc.errorcube,
                          header=spc.header,
                          normalize=False)

    cube_KZ_kwargs.update(lnKZ_xy_kwargs)

    hdu_K = cube_K(writeto=file_Ks, **cube_KZ_kwargs)
    hdu_Z = cube_Z(writeto=file_Zs, **cube_KZ_kwargs)

    return hdu_K, hdu_Z


if '__name__' == '__main__':
    make_KZ_cubes([0, 1, 2])
