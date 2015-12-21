"""
Match the SpIES catalog to the Tractor sweep files.
J. Moustakas - Siena College - 2015 May 21
"""

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from desispec.io.util import write_bintable
from astrometry.libkd.spherematch import match_radec as spherematch

cat = fits.getdata('cutout_350_365.fits',1)
sweepslice = ['000','001','002','003','004','005','350','351','352',
              '353','354','355','356','357','358','359']
sweeppath = os.path.join(os.getenv('DECALS_DIR'),'dr1','sweep')
sweepfile = [sweeppath+'/tractor-sweep-'+ss+'.fits' for ss in sweepslice]

tractor = None
for sw in sweepfile:
    print('Working on sweep {}'.format(sw))
    sweep = fits.getdata(sw,1)
    m1, m2, d12 = spherematch(sweep['ra'],sweep['dec'],
                              cat['ra'],cat['dec'],1.0/3600.0)
    print('Found {} matching objects'.format(len(m1)))
    if tractor is None:
        tractor = sweep[:][m1]
        spies = cat[:][m2]
    else: 
        tractor = vstack([Table(tractor),Table(sweep[:][m1])],join_type='exact')
        spies = vstack([Table(spies),Table(cat[:][m2])],join_type='exact')

spies.write('spies.fits')
tractor.write('tractor-spies.fits')

