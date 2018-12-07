""" Gets a quick fits file with "filtered" velocity differences """

import os
from numpy import nan
from astropy.io import fits

from config import chain_dir, file_mle_x2, file_Ks

Ks = fits.getdata(file_Ks)

mle2 = fits.getdata(file_mle_x2)
diff = mle2[10] - mle2[4]

# assure that we ran things in velocity separation space
assert diff[diff < 0].size == 0

diff[Ks[0] < 5] = nan
diff[Ks[1] < 5] = nan

hdu_diff = fits.PrimaryHDU(diff, header=fits.getheader(file_mle_x2))
hdu_diff.writeto(os.path.join(chain_dir, 'mle2_xoff_diff.fits'),
                 overwrite=True)

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
