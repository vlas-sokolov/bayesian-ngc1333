About
-----

This repository provides the code necessary to reproduce the
nested sampling results on the Green Bank Ammonia Survey ([GAS](https://arxiv.org/abs/1704.06318))
data on NGC 1333. It includes the scripts for data import, spectral
sampling, methods to parallelize the sampling on a spectral cube, and
scripts for the collection and visualization of the results.

Requirements
------------

The `requirements.txt` file contains a complete list
(`pip install -r requirements.txt` should do the job)

* [pyspecnest](https://github.com/vlas-sokolov/pyspecnest) is needed for sampling and cube-oriented results collection.
* [pymultinest](https://johannesbuchner.github.io/PyMultiNest/install) and [pyspeckit](https://github.com/pyspeckit/pyspeckit) as base dependencies (`pyspecnest` is essentially a wrapper of the two)
* Optional plotting dependencies used in imaging scripts:
    * [APLpy](https://github.com/aplpy/aplpy) for image plotting. A version of `1.1.1` (slightly older is fine too) is needed.
    * [ChainConsumer](https://github.com/Samreay/ChainConsumer) was used for generating corner plots.
    * [scikit-image](https://github.com/scikit-image/scikit-image) was used for small feature removal in plotting maps.

How it works
------------

The script hierarchy is as follows:

* A [configuration script](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/config.template.py) example should be copied (and adapted to your installation) into `config.py`. All folder paths, file names, tunable spectral model parameters and priors reside in it.
* A [cube loader script](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/opencube.py) is provided as well. Should be adapted to load a `pyspeckit`'s [`Cube`](https://pyspeckit.readthedocs.io/en/latest/cubes.html) or `CubeStack` instance via `opencube.make_cube` function.
    * Because loading the full FILE file into RAM is time-consuming, a `opencube.save_datacube` function has to be executed first, to prepare the spectra for lazy loading.
    * `import opencube; spc=opencube.make_cube(); opencube.save_datacube( spc)`
* A [script](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/innocent_script.py) preforming a sampling on a single spectrum can be invoked from command line with the following arguments:
    * `npeaks`: number of independent velocity components to sample with (defaults to `default_npeaks` in `config.py`)
    * `y`, `x`: Y/X pixel coordinates of a target spectrum (defaults to `default_yx` in `config.py`)
    * `plotting`: a Boolean specifying whether perform automatic plotting after the sampling has finished (defaults to 0: no plotting)
* A pool of processes is generated in the scheduling / parallelization [script](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/pool_xy.py), which can invoke a batch sampling run, for example, as follows:
    * `python pool_xy $NPEAKS $METHOD $CUT $NCPU` with the following arguments:
        * `$NPEAKS`: number of velocity components, defaults to 1.
        * `$METHOD`: `snr` or `Bfactor`, sets ordering of the process schedule (e.g., with `snr` high SNR regions are processed first), defaults to `snr`.
        * `$CUT`: stopping condition for the chosen method (e.g., method='snr' with cut=5 finishes when all SNR>5 pixels were sampled), defaults to 8.
        * `$N_CPU`: number of CPUs to use, defaults to 7.
    * `make_parcubes.py` and `make_KZ_cubes.py` are generating FITS files with MLE parameter estimates and maps of `lnK` and`lnZ`, respectively.
    * A sample plotting scripts of note is `make_presentable_Ks.py` (maps of Bayes factors, spectral multiplicity map), and `plot_spectra.py` (plots MLE parameters over selected spectra).

A note on parallelization: the sampling script is fully compatible with
Open MPI, meaning that one can run seven sampling processes in parallel
 by simply invoking:

`mpirun -np 7 python innocent_script.py $NPEAKS $Y $X`,

where `$NPEAKS` is a number of components to sample and `$X` and `$Y` are
pixel coordinates of the sampled spectrum. However, because all spectra are
treated independently, the sampling of the cube is an [embarrassingly parallel
problem](https://en.wikipedia.org/wiki/Embarrassingly_parallel), and the
standard Python `multiprocessing` library is used for thread spawning.

This code is released under an MIT open-source license. Contributions fitting the scope of the project are welcome, and bug reports can be filed in the Issues tab.
