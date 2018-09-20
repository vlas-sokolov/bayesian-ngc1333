from opencube import make_cube_shh
from pyspeckit.spectrum.units import SpectroscopicAxis
import numpy as np
from astropy import log
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
plt.rc('xtick', direction='in')
plt.rc('ytick', direction='in')
plt.rc('text', usetex=True)
plt.rc('font', **{'family' : "sans-serif"})
params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)

# the guys below are used to annotate npeaks map in make_presentable_Ks.py
from spectra_xy_list import xlist, ylist, labels

# will either make two subfigures (to be put on the sides of M-map)
# or makes a joint figure for spectra only for `split = False`
split = True

# colors for spectral components
cbrews = ['#e41a1c', '#377eb8', '#984ea3', '#4daf4a'][::-1]

spc = make_cube_shh()
spc.update_model("cold_ammonia")

mle_x1_file = 'nested-sampling/ngc1333-gas-mle-x1.fits'
mle_x2_file = 'nested-sampling/ngc1333-gas-mle-x2.fits'
Kfile = 'nested-sampling/ngc1333-Ks.fits'
Kcut = 5 # the heuristical ln(Z1/Z2) cut for model selection

x1_mle = fits.getdata(mle_x1_file)
x2_mle = fits.getdata(mle_x2_file)
Ks = fits.getdata(Kfile)

parcube = np.full_like(x2_mle, np.nan)
parcube[:6, :, :] = x1_mle[:6, :, :]
for i in range(1, 12):
    parcube[i, :, :][Ks[1]>5] = x2_mle[i, :, :][Ks[1]>5]

gridspec_kw={}
do_ave = False

if not split:
    plt.rc('font', size=10) # cosmetic tinkering
    fig, axarr = plt.subplots(xlist.size, 2, sharex=False,
                              sharey=False, figsize=(11.5, 4.6),
                              gridspec_kw=gridspec_kw)
else:
    plt.rc('font', size=14) # cosmetic tinkering
    f11, axarr11 = plt.subplots(xlist.size, 1, sharex=True,
                                sharey=False, figsize=(7.75, 7.6),
                                gridspec_kw=gridspec_kw)
    f22, axarr22 = plt.subplots(xlist.size, 1, sharex=True,
                                sharey=False, figsize=(3, 7.6),
                                gridspec_kw=gridspec_kw)
    axarr = np.array([axarr11, axarr22]).T

# make the high-res xarr for models
def highres_xarr(xarr, N, refX):
    # make sure the 11 and 22 cubes are in GHz
    xarr.convert_to_unit(u.GHz)
    xnparr = np.linspace(xarr.min(), xarr.max(), N)
    xarr_hires = SpectroscopicAxis(xarr=xnparr, refX=refX,
                    velocity_convention="radio")
    return xarr_hires

N = 1000 # resolution for modelled spectra
from pyspeckit.spectrum.models.ammonia_constants import freq_dict
spc.cubelist[0].xarr.velocity_convention = "radio"
spc.cubelist[1].xarr.velocity_convention = "radio"
xarr11, xarr22 = spc.cubelist[0].xarr, spc.cubelist[1].xarr
xarr11.refX = freq_dict["oneone"]*u.Hz
xarr22.refX = freq_dict["twotwo"]*u.Hz
h_xarr11 = highres_xarr(xarr11, N, freq_dict["oneone"]*u.Hz)
h_xarr22 = highres_xarr(xarr22, N, freq_dict["twotwo"]*u.Hz)
h_xarr11.convert_to_unit("km/s")
h_xarr22.convert_to_unit("km/s")

lab_fmt = r"$\sigma={}~\mathrm{{km~s^{{-1}}}}$"

for x, y, (ax11, ax22), lab in zip(xlist, ylist, axarr, labels):
    sp = spc.get_spectrum(x, y)
    sp11 = spc.cubelist[0].get_spectrum(x, y)
    sp22 = spc.cubelist[1].get_spectrum(x, y)
    sp11.xarr.convert_to_unit("km/s")
    sp22.xarr.convert_to_unit("km/s")
    pars = parcube[:, y, x]

    strK10 = r'$\ln K^{{\mathrm{{1}}}}_{{\mathrm{{0}}}} = {:.0f}$'.format(Ks[0, y, x])
    strK21 = r'$\ln K^{{\mathrm{{2}}}}_{{\mathrm{{1}}}} = {:.0f}$'.format(Ks[1, y, x])
    # plot the loc's of the spectra (cf. subsonic_analyses.py script)
    ax22.annotate(int(lab), xy=(0.85, 0.8), xycoords='axes fraction',
                  size='medium', bbox=dict(boxstyle="circle", fc="0.9"))
    ax11.text(x=0.85, y=0.80, s=strK10, transform=ax11.transAxes)
    ax11.text(x=0.85, y=0.58, s=strK21, transform=ax11.transAxes)
    ax11.text(x=0.5, y=0.80, s=r'$\mathrm{{(x, y) = ({}, {})}}$'.format(x, y),
              transform=ax11.transAxes, size='small')

    ax11.plot(sp11.xarr.value, sp11.data, drawstyle="steps-mid",
              color='black', lw=0.5, zorder=0.1)
    ax22.plot(sp22.xarr.value, sp22.data, drawstyle="steps-mid",
              color='black', lw=0.5, zorder=0.1)

    if np.isfinite(pars[6:18]).any():
        log.warn("Spectrum at (x,y)=({},{}) contains"
                 " an additional velocity component!".format(x, y))

    h_msum_11 = 0
    mcomps, h_11s, h_22s = [], [], []
    for p_ext, col in zip([pars[6:12], pars[0:6]], cbrews[2:]):
        try:
            h_11 = sp.specfit.get_model(h_xarr11, pars=p_ext)
            h_22 = sp.specfit.get_model(h_xarr22, pars=p_ext)
            # plot the individual LoS components
            ax11.plot(h_xarr11.value, h_11, lw=1.2, zorder=1, color=col)
            ax22.plot(h_xarr22.value, h_22, lw=1.2, zorder=1, color=col)
            h_11s.append(h_11)
            h_22s.append(h_22)
        except ValueError:
            pass

    # plot the total fit
    ax11.plot(h_xarr11.value, np.sum(h_11s, axis=0), lw=1.3,
              color='k', ls=':', zorder=1)
    ax22.plot(h_xarr22.value, np.sum(h_22s, axis=0), lw=1.3,
              color='k', ls=':', zorder=1)

    ax11.set_xlim(-15, 35)
    ax22.set_xlim(-2, 17)

    #ax22.set_ylim(-0.35, None)
    ax11.set_ylim(-0.3, None)
    if ax11.get_ylim()[1] < 0.5:
        ax11.set_ylim(None, 0.5)

    if not split and id(ax11) != id(axarr[-1, 0]):
        ax11.set_xticklabels([])
        ax22.set_xticklabels([])

    # don't just skip integer T_MB values
    ax11.yaxis.set_minor_locator(MultipleLocator(1))
    ax22.yaxis.set_minor_locator(MultipleLocator(0.5))
    # set the minor ticks equal to velocity reolution
    ax11.xaxis.set_minor_locator(MultipleLocator(1))
    ax22.xaxis.set_minor_locator(MultipleLocator(1))

ax11.set_xlabel(r"$\mathrm{Velocity~(km~s^{-1})}$")
ax22.set_xlabel(r"$\mathrm{Velocity~(km~s^{-1})}$")

if not split:
    top_ax11, top_ax22 = axarr[0]
else:
    top_ax11, top_ax22 = axarr11[0], axarr22[0]

top_ax11.text(0.04, 0.8, r'$\mathrm{NH_3~(1,1)}$',# size='large',
              transform=top_ax11.transAxes)
top_ax22.text(0.04, 0.8, r'$\mathrm{NH_3~(2,2)}$',# size='large',
              transform=top_ax22.transAxes)

# top left only
ax11.set_ylabel(r"$T_{\mathrm{MB}}~\mathrm{(K)}$")

if not split:
    fig.tight_layout()
    fig.savefig("figs/spectra.pdf", dpi=160)
else:
    f11.tight_layout()
    f22.tight_layout()
    f11.savefig("figs/spectra-11.pdf", dpi=160)
    f22.savefig("figs/spectra-22.pdf", dpi=160)

plt.show()
