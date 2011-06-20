pro photoztemplates_kcorrect, clobber=clobber
; jm10aug30ucsd - compute K-corrections with and without the zeropoint
; corrections 

    path = sings_path(/proj)+'photoztemplates/'
    phot = mrdfits(path+'phot.fits.gz',1)
    ngal = n_elements(phot)

    filters = photoztemplates_filterlist()
    nfilt = n_elements(filters)

    use = rebin(reform((k_lambda_eff(filterlist=filters) lt 5D4),nfilt,1),nfilt,ngal) ; <5 microns

    kcorr_emlines = photoztemplates_do_kcorrect(phot.z,phot.maggies,$
      phot.ivarmaggies*use,filterlist=filters,/emlines)
    im_mwrfits, kcorr_emlines, 'kcorr_emlines.fits', /clobber

    kcorr = photoztemplates_do_kcorrect(phot.z,phot.maggies,$
      phot.ivarmaggies*use,filterlist=filters)
    im_mwrfits, kcorr, 'kcorr.fits', /clobber

return
end
