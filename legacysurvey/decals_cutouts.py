from astropy.io import fits
import numpy as np
from desitarget import desi_mask

tt = fits.getdata('targets-dr2-0.2.2.fits')
bgs = np.where( (tt.DESI_TARGET & desi_mask.BGS_ANY) != 0 )[0]
url = "http://legacysurvey.org/viewer/?ra={ra}&dec={dec}&zoom=15&layer=decals-dr2p"

for i in bgs[0:100]:
    print url.format(ra=tt.RA[i], dec=tt.DEC[i])


    module load desitarget/0.2.0
    cd /project/projectdirs/desi/target
    select_targets /project/projectdirs/desiproc/dr2/sweep/ targets-dr2-0.2.0.fits -v

