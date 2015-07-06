#!/usr/bin/env python

def main():

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from astropy.io import fits
    from astropy.table import table
    import scipy.interpolate

    npix = 100
    
    edrpath = os.getenv('DECALS_DIR')

    hdu = fits.open('info.fits')
    #hdu = fits.open(os.path.join(edrpath,'info.fits'))
    info = hdu[1].data
    hdu.close()

    ra = info['RA']
    dec = info['DEC']
    dm_z = info['DM_Z']
    good = np.where(dm_z>-90)
    nobj = len(good)
    #print ra[good], dec[good], dm_z[good]

    xi = np.linspace(ra[good].min(),ra[good].max(),npix),
    yi = np.linspace(dec[good].min(),dec[good].max(),npix)
    #xi, yi = np.meshgrid(xi,yi)

    #Interpolate; there's also method='cubic' for 2-D data such as here
    zi = scipy.interpolate.griddata((ra[good],dec[good]),
            dm_z[good],(xi,yi),method='cubic',fill_value=0.0)
    print zi

    plt.figure(1)
    vmin = dm_z[good].min()
    vmax = dm_z[good].max()
    print vmin, vmax
    plt.imshow(zi,vmin=vmin,vmax=vmax,origin='lower',
            extent=[ra[good].min(),ra[good].max(),
            dec[good].min(),dec[good].max()],aspect='auto')
    cb = plt.colorbar()
    cb.set_clim(vmin=vmin,vmax=vmax)
    plt.savefig('test.png')
    plt.close()

if __name__ == '__main__':
    main()
