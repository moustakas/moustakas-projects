pro photoztemplates_write_kcorrect_models, out
; jm11may16ucsd - write out the models for M. Brown    

    vname = 'default.nolines'
    
    path = sings_path(/proj)+'photoztemplates/'
    phot = mrdfits(path+'phot.fits.gz',1)
    kcorr = mrdfits(path+'kcorr.fits.gz',1)
    ngal = n_elements(phot)
    
    k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
      lfile=lfile, vpath=vpath, vname=vname
    wave = k_lambda_to_centers(lambda)
    
    out = {galaxy: phot.galaxy, z: phot.z, wave: wave, $
      flux: vmatrix#kcorr.k_coeffs}
    im_mwrfits, out, path+'kcorrect_models.fits', /clobber
    
return
end
    
