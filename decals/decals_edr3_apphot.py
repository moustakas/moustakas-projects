if __name__ == '__main__':
"""
Perform aperture photometry on EDR3 images and residual images.
J. Moustakas, 2013 Mar 05, Siena
"""

    import os
    import photutils
    import numpy as np
    from glob import glob
    from astropy.io import fits
    from astropy.table import table

    edrpath = '/Users/ioannis/research/projects/decals/edr3/'
    #edrpath = '/global/work/decam/cats/edr3/'

    band = ['g','r','z']

    maxrad = 15.0 # maximum aperture radius [pixels]
    radius = np.arange(0.0,maxrad)+1.0 # [pixels]

    allbrick = glob(edrpath+'tractor/*')
    nbrick = len(allbrick)
    for ibrick in allbrick:
        bigbrick = os.path.basename(ibrick)
        catfile = glob(os.path.join(edrpath,'tractor',bigbrick)+'/tractor-*.fits')
        ncat = len(catfile)

        for icat in catfile:
            outcatfile = icat.replace('tractor','tractor2')
            hducat = fits.open(icat)
            cat = hducat[1].data
            hducat.close()
            nobj = len(cat)

            xypos = np.array([cat['BX'],cat['BY']]).transpose() # pixel positions

            brickname = cat[0]['brickname']
            coaddpath = os.path.join(edrpath,'coadd',bigbrick,brickname)
            for iband in band:
                hduimage = fits.open(coaddpath+'/decals-'+brickname+'-image-'+iband+'.fits')
                hduinvvar = fits.open(coaddpath+'/decals-'+brickname+'-invvar-'+iband+'.fits')
                hdumodel = fits.open(coaddpath+'/decals-'+brickname+'-model-'+iband+'.fits')

                image = hduimage[0].data
                invvar = hduinvvar[0].data
                model = hdumodel[0].data
                resid = image - model # residual image

                with np.errstate(divide='ignore'):
                    imsigma = 1.0/np.sqrt(invvar)
                    imsigma[invvar == 0] = 0

                # get the aperture photometry on the data and residual
                # images
                flux = []
                flux_resid = []
                for rad in radius:
                    aper = photutils.CircularAperture(xypos,rad)
                    flux.append(photutils.aperture_photometry(image,aper,error=imsigma))
                    flux_resid.append(photutils.aperture_photometry(resid,aper,error=imsigma))

                # stack the table of fluxes
                apphot = hstack(flux)
                apphot_resid = hstack(flux_resid)

                # combine the aperture photometry with CAT and write
                # out! not sure how to do this yet!

                # clean up
                hduimage.close()
                hduinvvar.close()
                hdumodel.close()
