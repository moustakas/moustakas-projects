pro sdss_lowz_isedfit, models=models, isedfit=isedfit, clobber=clobber
; jm11mar24ucsd - 

    iopath = sdss_path()+'lowz_isedfit/'
    paramfile = iopath+'sdss_lowz_isedfit.par'

; ---------------------------------------------------------------------------
; build the model grids    
    if keyword_set(models) then isedfit_models, paramfile, $
      iopath=iopath, clobber=clobber

; ---------------------------------------------------------------------------
; do the SED fitting
    if keyword_set(isedfit) then begin
       lowz = mrdfits(sdss_path(/lowz)+'lowz_catalog.dr6.fits.gz',1)
       sdss_to_maggies, maggies, ivarmaggies, calibobj=lowz, flux='model'

       isedfit, paramfile, maggies, ivarmaggies, lowz.zdist, result, $
         iopath=iopath, outprefix=outprefix, galchunksize=1000, $
         clobber=clobber, debug=debug, index=index
    endif

return
end
