from astropy.io import fits

mle_x2_file = 'nested-sampling/ngc1333-gas-mle-x2.fits'
Ks = fits.getdata('nested-sampling/ngc1333-Ks.fits')

mle2 = fits.getdata(mle_x2_file)
diff = mle2[10] - mle2[4]

# assure that we ran things in velocity separation space
assert diff[diff < 0].size == 0

from numpy import nan
diff[Ks[0] < 5] = nan
diff[Ks[1] < 5] = nan

hdu_diff = fits.PrimaryHDU(diff, header=fits.getheader(mle_x2_file))
hdu_diff.writeto('nested-sampling/mle2_xoff_diff.fits', overwrite=True)

#plt.figure()
#diff_good = diff.copy()
#diff_okayiguess = diff.copy()
#diff_bad = diff.copy()
#diff_good[(Ks[0] < 5) | (Ks[1] < 5)] = np.nan
#diff_okayiguess[(Ks[0] < 0) | (Ks[1] < 0)] = np.nan
#plt.imshow(diff_bad, origin='lower', vmax=2, alpha=0.1)
#plt.imshow(diff_okayiguess, origin='lower', vmax=2, alpha=0.6)
#plt.imshow(diff_good, origin='lower', vmax=2, alpha=1)
#plt.colorbar()
