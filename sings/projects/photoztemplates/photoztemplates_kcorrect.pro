pro photoztemplates_kcorrect, clobber=clobber
; jm10aug30ucsd - compute K-corrections with and without the zeropoint
; corrections 

    path = sings_path(/proj)+'photoztemplates/'
    phot = mrdfits(path+'phot.fits.gz',1)
    ngal = n_elements(phot)

    filters = photoztemplates_filterlist()
    nfilt = n_elements(filters)

    use = intarr(nfilt,ngal)+1
    use[12:16,*] = 0

    kcorr = photoztemplates_do_kcorrect(phot.z,phot.maggies,$
      phot.ivarmaggies*use,filterlist=filters)

    im_mwrfits, kcorr, 'kcorr.fits', /clobber

return
end
