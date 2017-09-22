#!/usr/bin/env python

"""Get the average seeing in the DEEP2 16- and 23-hour fields."""

# match, repstr(file_basename(strtrim(jj.cpimage,2)),'.fz','')+strtrim(jj.extname,2), strtrim(bb.filename,2)+strtrim(bb.ccdname,2), m1, m2

from __future__ import division, print_function

import os
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from astropy.io import fits

from projects.desi.common import *

# Global variables.
dr1dir = os.getenv('DECALS_DIR_DR1')
dr1mdir = dr1dir+'m'

def get_ccdinfo(brickwcs=None,deep2=None):
    """Get info on this brick and on the CCDs touching it."""

    if deep2 is not None:
        ccdfile = os.path.join(dr1mdir,'deep2-ccds-zeropoints.fits')
    else:
        ccdfile = os.path.join(dr1dir,'decals-ccds-zeropoints.fits')
    allccdinfo =  fits_table(ccdfile)

    # Get all the CCDs that touch this brick
    these = ccds_touching_wcs(brickwcs,allccdinfo)
    ccdinfo = allccdinfo[these]
    
    return ccdinfo

def main():

    pixscale = 1.0/3600.0 # [deg/pix]

    # Field 2
    ra = 252.45268
    dec = 34.942833
    dra = 1.75 # deg
    ddec = 0.5 # deg

    W = dra/pixscale
    H = ddec/pixscale
    wcs_field2 = Tan(ra, dec, W/2.+0.5, H/2.+0.5, 
                     -pixscale, 0., 0., pixscale,
                     float(W), float(H))
    ccds_field2 = get_ccdinfo(wcs_field2,deep2=1)

    # Field 3
    ra = 352.45486
    dec = 0.15857297
    dra = 2.25 # deg
    ddec = 0.55 # deg

    W = dra/pixscale
    H = ddec/pixscale
    wcs_field3 = Tan(ra, dec, W/2.+0.5, H/2.+0.5, 
                     -pixscale, 0., 0., pixscale,
                     float(W), float(H))
    ccds_field3 = get_ccdinfo(wcs_field3)

    bins = 25
    srange = [0.6,2.2]

    sns.set(style='white',font_scale=1.4)

    fig = plt.figure(figsize=(8,6))
    plt.hist(ccds_field2.seeing*0.262/0.27,bins,range=srange,
             normed=False,color="#6495ED",alpha=0.5,
             label='DEEP2/Field 2 (16h)')
    plt.hist(ccds_field3.seeing*0.262/0.27,bins,range=srange,
             normed=False,color="#F08080",alpha=0.4,
             label='DEEP2/Field 3 (23h)')
    plt.xlabel('FWHM Seeing (arcsec)')
    plt.ylabel('Number of CCDs')
    plt.xlim([0.6,2.1])
    plt.legend()
    fig.subplots_adjust(bottom=0.2)
    plt.savefig('decals_deep2_seeing.png')

if __name__ == "__main__":
    main()
