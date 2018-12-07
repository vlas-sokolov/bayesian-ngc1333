import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
from matplotlib.colors import ListedColormap
import aplpy

# sane defaults for aplpy-generated figures
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
plt.rc('text', usetex=True)
plt.rc('figure', figsize=(4, 6))
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='STIXGeneral')


def header_flatten(head):
    """ Flattens 3D header into 2D. Surely it exists somewhere already... """
    flathead = head.copy()
    for key in flathead.keys():
        if key.endswith('3'):
            flathead.pop(key)
    flathead['NAXIS'] = 2
    flathead['WCSAXES'] = 2

    return flathead


Kcut = 5 # the heuristical ln(Z1/Z2) cut for model selection

###########~~~~~~~~~~~ Load in GAS DR1 results on NGC1333 ~~~~~~~~~~###########
mle_x2_file = 'nested-sampling/ngc1333-gas-mle-x2.fits'
mle_x1_file = 'nested-sampling/ngc1333-gas-mle-x1.fits'
Ks = fits.getdata('nested-sampling/ngc1333-Ks.fits')

mle2 = fits.getdata(mle_x2_file)
mle1 = fits.getdata(mle_x1_file)

npeaks_map = np.full_like(Ks[0], np.nan)
xoff_mle = np.full_like(Ks[0], np.nan)

sig_mle = np.full_like(Ks[0], np.nan)
sig_mle[Ks[0] > Kcut] = mle1[3][Ks[0] > Kcut]

# results from pyspeckit (i.e. DR1 release)
sig_dr1 = fits.getdata('gasdata/NGC1333_Sigma_DR1_rebase3_flag.fits')
esig_dr1 = fits.getdata('gasdata/NGC1333_eSigma_DR1_rebase3_flag.fits')

sig_dr1[sig_dr1==0] = np.nan
###########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###########


# make the ln(K)>Kcut based map of LoS component numbers
npeaks_map[Ks[0] <= Kcut] = 0
npeaks_map[Ks[0] > Kcut] = 1
npeaks_map[(Ks[0] > Kcut) & (Ks[1] > Kcut)] = 2
hdu_npeaks = fits.PrimaryHDU(npeaks_map,
                header=header_flatten(fits.getheader(mle_x1_file)))
hdu_npeaks.writeto('nested-sampling/npeaks_cut5.fits', overwrite=True)
fig = aplpy.FITSFigure(hdu_npeaks)
cmaplst = ['#f1f7d2', '#7fcdbb', '#2c7fb8'] # stolen from colorbrewer2
lcmap = ListedColormap(cmaplst)
fig.show_colorscale(cmap=lcmap, vmin=-0.5, vmax=2.5)
fig._ax2.xaxis.set_tick_params('both', color='black')
fig._ax2.yaxis.set_tick_params('both', color='black')
fig._ax1.xaxis.set_tick_params('both', color='black')
fig._ax1.yaxis.set_tick_params('both', color='black')
fig.show_colorbar()
fig.colorbar.set_ticks([0, 1, 2])
fig.show_contour('gasdata/NGC1333_Sigma_DR1_rebase3_flag.fits', levels=[0],
                 linewidths=0.5, colors='black')
fig.savefig("npeaks_map.png", dpi=140)

# the ones below are the "overlapping" & "good" values
sig_mle_x1 = sig_mle[np.isfinite(sig_dr1) & (npeaks_map==1)]
sig_dr1_x1 = sig_dr1[np.isfinite(sig_dr1) & (npeaks_map==1)]
# the DR1 results that were found to have two LoS components
sig_dr1_x2_bad = sig_dr1[np.isfinite(sig_dr1) & (npeaks_map==2)]
# the MLE sigma values where two components are
sig_mle_x2_good = np.hstack([mle2[3][npeaks_map==2], mle2[3][npeaks_map==2]])
# the MLE sigma values where two components are
sig_mle_x1_good = mle1[3][npeaks_map==1]
sig_mle_all_good = np.hstack([sig_mle_x1_good, sig_mle_x2_good])
# all MLE results but depending on whether or not they were in GAS DR1
sig_mle_over_x1_good = mle1[3][np.isfinite(sig_dr1) & (npeaks_map==1)]
sig_mle_over_x2_good = np.hstack([mle2[3][np.isfinite(sig_dr1) & (npeaks_map==2)],
                                  mle2[3][np.isfinite(sig_dr1) & (npeaks_map==2)]])
sig_mle_overall_good = np.hstack([sig_mle_over_x1_good, sig_mle_over_x2_good])


# the overlapping values, where we're sure only one component is present,
# In [84]: pearsonr(sig_mle_x1, sig_dr1_x1)
# Out[84]: (0.9468249988181304, 0.0)
# ... are consistent between the two methods
# Would be interesting to see where the outliers are!

from matplotlib import pyplot as plt
histkwargs = dict(bins=80, alpha=0.5, range=(0, 1.4))

plt.figure("How are the DR1 results biased by 2nd components?",
           figsize=(11, 5))
plt.hist([sig_dr1_x1, sig_dr1_x2_bad], stacked=True,
         label='Values fit on 2-component spectra taken out', **histkwargs)
plt.legend()
plt.show()

# Okay, so what does this one show?
# Spectra with npeaks=1 but not npeaks=2 detection,
# that overlap on both GAS DR1 and my MLE results
from scipy.stats import pearsonr
plt.figure("Overlapping values shown: Pearson's"
           " r={:.4}".format(pearsonr(sig_dr1_x1, sig_mle_x1)[0]),
           figsize=(11, 5))
plt.hist(sig_dr1_x1, label='GAS DR1', **histkwargs)
plt.hist(sig_mle_x1, label='GAS MLE', **histkwargs)
plt.legend()
plt.text(1.02, 500, 'Edge of the prior, sorry!', rotation=90)
plt.axvline(1.0, ls=':', color='k')
plt.show()

plt.figure("All MLE results, for npeaks=1 and npeaks=2",
           figsize=(11, 5))
plt.hist([sig_mle_x1_good, sig_mle_x2_good], stacked=True,
         label='GAS MLE all', **histkwargs)
plt.show()


plt.figure("What about extra values fit in lower S/N regions?",
                   figsize=(11, 5))
plt.hist(sig_mle_all_good, label='All data', **histkwargs)
plt.hist(sig_mle_overall_good, label='Overlapping', **histkwargs)
plt.legend()
plt.show()
