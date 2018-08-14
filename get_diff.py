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
