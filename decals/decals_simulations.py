#!/usr/bin/env python

import numpy as np
import galsim

def makestamps(nstamp=50):
    """Document me!

    TODO: Include 2-component galaxies.
    """
    nmin = 0.8
    nmax = 5.0
    r50min = 0.5
    r50max = 2.0
    bamin = 0.1
    bamax = 1.0

    info = dict(
        SERSICN = np.random.uniform(nmin,nmax,nstamp),
        R50 = np.random.uniform(r50min,r50max,nstamp),
        BA = np.random.uniform(bamin,bamax,nstamp),
        PHI = np.random.uniform(0.0,360,nstamp)
        )

    for iobj in range(nstamp):
        gal = galsim.Sersic(info['SERSICN'][iobj],half_light_radius=info['R50'][iobj])
        gal = gal.shear(q=info['BA'][iobj],beta=info['PHI'][iobj]*galsim.degrees)
        im = gal.drawImage()
        galsim.fits.write(im,file_name='junk.fits')
        

def makeimages():

    help(galsim.Sersic)
           >>> im1 = bulge.drawImage()
               >>> im2 = disk.drawImage(image=im1, add_to_image=True)
               >>> assert im1 is im2
       
               >>> full_image = galsim.Image(2048, 2048, scale=pixel_scale)
               >>> b = galsim.BoundsI(x-32, x+32, y-32, y+32)
               >>> stamp = obj.drawImage(image = full_image[b])
               >>> assert (stamp.array == full_image[b].array).all()



def main():

    import os
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='DECaLS simulations.')

    parser.add_argument('--makestamps', action='store_false', 
                        help='make postage stamps')

    
    

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
