pro bcgs_sdss_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  chi2grid_qaplot=chi2grid_qaplot, measure=measure, clobber=clobber, $
  debug=debug
; jm10jul23 - measure stellar masses for the SDSS comparison sample of BCGs
; jm10sep28 - updated to a new SDSS sample
    
    iopath = ages_path(/projects)+'bcgs/'
    paramfile = iopath+'bcgs_sdss_isedfit.par'
    phot = mrdfits(iopath+'sdss_sample_v1.fits.gz',1)
;   phot = mrdfits(iopath+'sdss_08stott.fits.gz',1)

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       isedfit, paramfile, phot.maggies, phot.ivarmaggies, phot.z, $
         iopath=iopath, outprefix=outprefix, clobber=clobber, $
         debug=debug, index=index, write_chi2grid=1
    endif 

; --------------------------------------------------
; make some QAplots
    if keyword_set(qaplot) then begin
       isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
         index=index, clobber=clobber, outprefix=outprefix
    endif

    if keyword_set(chi2grid_qaplot) then begin
       isedfit_chi2grid_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
         index=index, clobber=clobber, outprefix=outprefix
    endif

; --------------------------------------------------
; measure rest-frame quantities
    if keyword_set(measure) then begin
       isedfit_measure, paramfile, measure, isedfit, iopath=iopath, $
         clobber=clobber, outprefix=outprefix
    endif

return
end
