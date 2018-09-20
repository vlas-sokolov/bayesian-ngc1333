""" Generating figures for the upcoming letter """
import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
from matplotlib.colors import ListedColormap
import aplpy
from skimage import morphology

# sane defaults for aplpy-generated figures
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
plt.rc('text', usetex=True)
plt.rc('figure', figsize=(4, 6))
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='STIXGeneral')

mle_x1_file = 'nested-sampling/ngc1333-gas-mle-x1.fits'
gascontourfile = 'gasdata/NGC1333_Sigma_DR1_rebase3_flag.fits'
Kfile = 'nested-sampling/ngc1333-Ks.fits'
Kcut = 5 # the heuristical ln(Z1/Z2) cut for model selection

def header_flatten(head):
    """ Flattens 3D header into 2D. Surely it exists somewhere already... """
    flathead = head.copy()
    for key in flathead.keys():
        if key.endswith('3'):
            flathead.pop(key)
    flathead['NAXIS'] = 2
    flathead['WCSAXES'] = 2

    return flathead

def show_filtered_contours(fig, Karr, Kcut=Kcut, skimage_kwargs={},
        header=header_flatten(fits.getheader(Kfile)), **kwargs):
    """ Overplots contours with small features removed """
    min_size = skimage_kwargs.pop('min_size', 3)
    levels = kwargs.pop('level', [0.5])
    Karr_clean = morphology.remove_small_objects(Karr > Kcut,
                        min_size=min_size).astype(int)
    hdu_clean = fits.PrimaryHDU(Karr_clean,
            header=header_flatten(fits.getheader(Kfile)))
    fig.show_contour(data=hdu_clean, levels=levels, **kwargs)

def tick_vaccination(fig):
    """ Fixes x and y ticks on aplpy FITSFigures in mpl 2.0+ """
    fig._ax2.yaxis.set_tick_params('both', color='black')
    fig._ax2.xaxis.set_tick_params('both', color='black')
    fig._ax1.yaxis.set_tick_params('both', color='black')
    fig._ax1.xaxis.set_tick_params('both', color='black')

Ks = fits.getdata(Kfile)

gas_contour_kwargs = dict(levels=[0], linewidths=0.5,
        colors='#ef6548', linestyles='--')
kcut_contour_kwargs = dict(linewidths=0.5, colors='white')

hdu_K01 = fits.PrimaryHDU(Ks[0], header=header_flatten(fits.getheader(Kfile)))
fig01 = aplpy.FITSFigure(hdu_K01)
fig01.show_colorscale(cmap='viridis', vmin=-2, stretch='power', exponent=0.3)
fig01.set_nan_color('#f0f0f0')
tick_vaccination(fig01)
fig01.show_colorbar()
fig01.show_contour(gascontourfile, **gas_contour_kwargs)
show_filtered_contours(fig01, Ks[0], **kcut_contour_kwargs)
fig01.colorbar.set_axis_label_text(r'$\ln K^{\mathrm{1}}_{\mathrm{0}}$')
fig01.savefig("figs/Ks_01_map.pdf", dpi=140)

hdu_K21 = fits.PrimaryHDU(Ks[1],
                header=header_flatten(fits.getheader(Kfile)))
fig21 = aplpy.FITSFigure(hdu_K21)
fig21.show_colorscale(cmap='viridis', vmin=-2, stretch='power', exponent=0.5)
fig21.set_nan_color('#f0f0f0')
tick_vaccination(fig21)
fig21.show_colorbar()
fig21.show_contour(gascontourfile, **gas_contour_kwargs)
show_filtered_contours(fig21, Ks[1], **kcut_contour_kwargs)
fig21.colorbar.set_axis_label_text(r'$\ln K^{\mathrm{2}}_{\mathrm{1}}$')
fig21.savefig("figs/Ks_21_map.pdf", dpi=140)

# not sure which one is better... so why not both?
for flag_noise, suffix in zip([True, False], ['_clean', '_all']):
    mle1 = fits.getdata(mle_x1_file)
    min_size = 9
    if flag_noise:
        Karr_clean = morphology.remove_small_objects(Ks[0] > 0,
                        min_size=min_size)
        mle1[4][(Ks[0] < 0) | ~Karr_clean] = np.nan

    hdu_xoff = fits.PrimaryHDU(mle1[4],
                    header=header_flatten(fits.getheader(Kfile)))
    figxoff = aplpy.FITSFigure(hdu_xoff)
    figxoff.show_colorscale(cmap='RdYlBu_r')
    figxoff.set_nan_color('#f0f0f0')
    tick_vaccination(figxoff)
    figxoff.show_colorbar()
    gas_contour_kwargs.update({'linewidths': 0.7}) # hey that's not cool!
    figxoff.show_contour(gascontourfile, **gas_contour_kwargs)
    show_filtered_contours(figxoff, Ks[0], Kcut=0,
            skimage_kwargs=dict(min_size=min_size),
            linewidths=0.7, colors='black')
    figxoff.colorbar.set_axis_label_text(
            r'$V_{\mathrm{lsr}}~\mathrm{(km~s^{-1})}$')
    figxoff.savefig("figs/mle_xoff{}.pdf".format(suffix), dpi=140)

# make the ln(K)>Kcut based map of LoS component numbers
npeaks_map = np.full_like(Ks[0], np.nan)
npeaks_map[Ks[0] <= Kcut] = 0
npeaks_map[Ks[0] > Kcut] = 1
Karr_clean = morphology.remove_small_objects(Ks[0] > Kcut,
        min_size=3).astype(int)
npeaks_map[(Karr_clean==0) & np.isfinite(npeaks_map)] = 0
npeaks_map[(Ks[0] > Kcut) & (Ks[1] > Kcut)] = 2
hdu_npeaks = fits.PrimaryHDU(npeaks_map,
                header=header_flatten(fits.getheader(mle_x1_file)))
hdu_npeaks.writeto('nested-sampling/npeaks_cut5.fits', overwrite=True)
fig = aplpy.FITSFigure(hdu_npeaks)
cmaplst = ['#f1f7d2', '#7fcdbb', '#2c7fb8'] # stolen from colorbrewer2
lcmap = ListedColormap(cmaplst)
fig.show_colorscale(cmap=lcmap, vmin=-0.5, vmax=2.5)
tick_vaccination(fig)
fig.show_colorbar()
fig.colorbar.set_ticks([0, 1, 2])
fig.show_contour(gascontourfile, **gas_contour_kwargs)

# nothing fancy follows - we've got to fine-tune all the arrows manually!
from spectra_xy_list import xlist, ylist, labels
ax = fig._ax1
x1, x2, x3, x4, x5 = xlist
y1, y2, y3, y4, y5 = ylist
l1, l2, l3, l4, l5 = labels

ax.annotate(l1, xy=(x1+1, y1+1), xycoords='data',
            xytext=(+30, -10), textcoords='offset points',
            bbox=dict(boxstyle="circle", fc="0.8", alpha=0.8),
            arrowprops=dict(arrowstyle="-", shrinkA=0, shrinkB=0,
                            connectionstyle="arc3, rad=0.3",
                            alpha=0.8))
ax.annotate(l2, xy=(x2+1, y2+1), xycoords='data',
            xytext=(30, -25), textcoords='offset points',
            bbox=dict(boxstyle="circle", fc="0.8", alpha=0.8),
            arrowprops=dict(arrowstyle="-", shrinkA=0, shrinkB=0,
                            connectionstyle="arc3, rad=-0.4",
                            alpha=0.8))
ax.annotate(l3, xy=(x3+1, y3+1), xycoords='data',
            xytext=(-60, -40), textcoords='offset points',
            bbox=dict(boxstyle="circle", fc="0.8", alpha=0.8),
            arrowprops=dict(arrowstyle="-", shrinkA=0, shrinkB=0,
                            connectionstyle="arc3, rad=0.4",
                            alpha=0.8))
ax.annotate(l4, xy=(x4+1, y4+1), xycoords='data',
            xytext=(+60, 40), textcoords='offset points',
            bbox=dict(boxstyle="circle", fc="0.8", alpha=0.8),
            arrowprops=dict(arrowstyle="-", shrinkA=0, shrinkB=0,
                            connectionstyle="arc3, rad=+0.2",
                            alpha=0.8))
ax.annotate(l5, xy=(x5+1, y5+1), xycoords='data',
            xytext=(-60, 40), textcoords='offset points',
            bbox=dict(boxstyle="circle", fc="0.8", alpha=0.8),
            arrowprops=dict(arrowstyle="-", shrinkA=0, shrinkB=0,
                            connectionstyle="arc3, rad=-0.2",
                            alpha=0.8))



from matplotlib.patches import Circle
for (x, y) in zip(xlist, ylist):
    r = Circle((x+1, y+1), radius=.5, zorder=1, fc="0.8", ec="k", lw=0.2)
    ax.add_patch(r)
fig.savefig("figs/npeaks_map.pdf", dpi=140)

plt.show()
