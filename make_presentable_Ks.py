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


Kfile = 'nested-sampling/ngc1333-Ks.fits'

Kcut = 5 # the heuristical ln(Z1/Z2) cut for model selection

Ks = fits.getdata(Kfile)

hdu_K01 = fits.PrimaryHDU(Ks[0], header=header_flatten(fits.getheader(Kfile)))
fig01 = aplpy.FITSFigure(hdu_K01)
fig01.show_colorscale(cmap='viridis', vmin=-2, stretch='power', exponent=0.3)
fig01._ax2.yaxis.set_tick_params('both', color='black')
fig01._ax2.yaxis.set_tick_params('both', color='black')
fig01._ax1.yaxis.set_tick_params('both', color='black')
fig01._ax1.yaxis.set_tick_params('both', color='black')
fig01.show_colorbar()
fig01.show_contour('gasdata/NGC1333_Sigma_DR1_rebase3_flag.fits', levels=[0],
                   linewidths=0.5, colors='white')
fig01.show_contour(hdu_K01, levels=[5], linewidths=0.5, colors='red')
fig01.colorbar.set_axis_label_text(r'$\ln K_{\mathrm{10}}$')
fig01.savefig("Ks_01_map.png", dpi=140)

hdu_K21 = fits.PrimaryHDU(Ks[1],
                header=header_flatten(fits.getheader(Kfile)))
fig21 = aplpy.FITSFigure(hdu_K21)
fig21.show_colorscale(cmap='viridis', vmin=-2, stretch='power', exponent=0.5)
fig21._ax2.yaxis.set_tick_params('both', color='black')
fig21._ax2.yaxis.set_tick_params('both', color='black')
fig21._ax1.yaxis.set_tick_params('both', color='black')
fig21._ax1.yaxis.set_tick_params('both', color='black')
fig21.show_colorbar()
fig21.show_contour('gasdata/NGC1333_Sigma_DR1_rebase3_flag.fits', levels=[0],
                 linewidths=0.5, colors='white')
fig21.show_contour(hdu_K21, levels=[5], linewidths=0.5, colors='red')
fig21.colorbar.set_axis_label_text(r'$\ln K_{\mathrm{21}}$')
fig21.savefig("Ks_21_map.png", dpi=140)
