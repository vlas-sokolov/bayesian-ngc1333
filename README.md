About
-----

This repository provides the code necessary to reproduce the
nested sampling results on the Green Bank Ammonia Survey ([GAS](https://arxiv.org/abs/1704.06318))
data on NGC 1333 (Sokolov et al. 2020, [ApJL](https://iopscience.iop.org/article/10.3847/2041-8213/ab8018) / [arXiv](https://arxiv.org/abs/2003.07644)). It includes the scripts for data import, spectral
sampling, methods to parallelize the sampling on a spectral cube, and
scripts for the collection and visualization of the results.

Evidence / Bayes factor maps, as well as MLEs of parameter estimation are available at https://doi.org/10.7910/DVN/PDH6VT.

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

Usage instructions
------------------

To avoid playing a lengthy game of whack-an-error, please read the following items carefully. The following describes the script hierarchy and a rough order of things to do before the sampler is ready to be used:

1. Install the requirements above. Make sure you can run MultiNest via PyMultiNest - try an example [here](https://johannesbuchner.github.io/PyMultiNest/install.html#running-some-code) to make sure everything works.
2. Set up a configuration file by copying a [template file](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/config.template.py), which should be adapted for your local setup and named as `config.py`. All folder paths, file names, tunable spectral model parameters and priors reside in it.
3. Pull the GAS data files, both the [FITS cubes](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/IP1ZZL) as well as velocity [dispersion files](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/YAXXAY) (for paper figures only), and place them into the `data/` folder.
4. Digest the GAS data for lazy loading of the spectra. A [cube loader script](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/opencube.py) is provided for it. It should be adapted to load a `pyspeckit`'s [`Cube`](https://pyspeckit.readthedocs.io/en/latest/cubes.html) or `CubeStack` instance via `opencube.make_cube` function and ran as follows:
    * Because loading the full FILE file into RAM is time-consuming, a `opencube.save_datacube` function has to be executed first, to prepare the spectra for lazy loading.
    * Make sure `nested-sampling/cubexarr/` folder exists (or whichever folder is defined as `cube_storage_dir` variable in `config.py`).
    * Run `import opencube; spc=opencube.make_cube(); opencube.save_datacube(spc)`.
5. A [script](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/innocent_script.py) that runs a **sampling on a single spectrum** can be invoked from command line with the following arguments:
    * `npeaks`: number of independent velocity components to sample with (defaults to `default_npeaks` in `config.py`)
    * `y`, `x`: Y/X pixel coordinates of a target spectrum (defaults to `default_yx` in `config.py`)
    * `plotting`: a Boolean specifying whether to perform automatic plotting after the sampling has finished (defaults to 0: no plotting)
6. A pool of processes is generated in the scheduling / parallelization [script](https://github.com/vlas-sokolov/bayesian-ngc1333/blob/master/pool_xy.py), which can invoke a **batch sampling on a full spectral cube**, for example, as follows:
    * `python pool_xy.py $NPEAKS $METHOD $CUT $NCPU` with the following arguments:
        * `$NPEAKS`: number of velocity components, defaults to 1.
        * `$METHOD`: `snr` or `Bfactor`, sets ordering of the process schedule (e.g., with `snr` high SNR regions are processed first), defaults to `snr`.
        * `$CUT`: stopping condition for the chosen method (e.g., method='snr' with cut=5 finishes when all SNR>5 pixels were sampled), defaults to 8.
        * `$N_CPU`: number of CPUs to use, defaults to 7.
7. Finally, to make evidence and Bayes factor FITS maps, run `make_KZ_cubes.py`; and to make "best-fit", i.e., MLE point estimate parameter cubes, run `make_parcubes.py`.
8. Post-processing and fancy plotting scripts are inside `make_presentable_Ks.py` (maps of Bayes factors, spectral multiplicity map), and `plot_spectra.py` (plots MLE parameters over selected spectra) files.

A note on parallelization
-------------------------

The sampling script is fully compatible with
Open MPI, meaning that one can run seven sampling processes in parallel
 by simply invoking:

`mpirun -np 7 python innocent_script.py $NPEAKS $Y $X`,

where `$NPEAKS` is a number of components to sample and `$X` and `$Y` are
pixel coordinates of the sampled spectrum. However, because all spectra are
treated independently, the sampling of the cube is an [embarrassingly parallel
problem](https://en.wikipedia.org/wiki/Embarrassingly_parallel), and the
standard Python `multiprocessing` library is used for thread spawning in `pool_xy.py`.

FAQ
---

> `FileNotFoundError: [Errno 2] No such file or directory: '/path-to-project-folder/nested-sampling/cubexarr/ngc1333-gas-xarr.npy'`

Make sure you follow the instructions on point (4) above. The scripts are not smart enough (feel free to make a PR for that) to do it automatically on the first run.

> `Fortran runtime error: File already opened in another unit` or other errors of the can't-find-chains-files type

Unfortunately, it seems that the current MutliNest version (`LEN=100` in src files [here](https://github.com/JohannesBuchner/MultiNest/blob/master/src/nested.F90)) has a hard limit of 100 characters on its *absolute path* variables, and will truncate it to a hundred chars in a criminally silent manner. I found it out the hard way, so make sure you don't run over it.

> `too many files open` error

Routinely happens on some platforms. Follow instructions [here](https://superuser.com/questions/1200539/cannot-increase-open-file-limit-past-4096-ubuntu) or [here](https://stackoverflow.com/questions/34588/how-do-i-change-the-number-of-open-files-limit-in-linux) to increase the maximum number of files your OS lets you open at the same time.

*This code is released under an MIT open-source license. Contributions fitting the scope of the project are welcome, and bug reports can be filed in the Issues tab.*
