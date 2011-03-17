pro bcgs_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  chi2grid_qaplot=chi2grid_qaplot, measure=measure, clobber=clobber, $
  debug=debug
; jm10jun02 - measure stellar masses for A. Gonzalez's sample of BCGs

    iopath = ages_path(/projects)+'bcgs.test/'
    paramfile = iopath+'bcgs_isedfit.par'

; read the photometry and convert to maggies; do not use ch4 when it
; shifts redward of 3 microns, and do not use the Bw-band blueward of
; ~2200 A
    vv = 'v3'
    phot = mrdfits(iopath+'bcgs_photometry_'+vv+'.fits.gz',1)
    sample = rsex(iopath+'bcgs_sample_'+vv+'.sex')
    redshift = sample.z
    bootes_to_maggies, phot, maggies, ivarmaggies, $
      use_aper='08', filterlist=filterlist

    keep = where((strmatch(filterlist,'*ufilter*') eq 0) and $
      (strmatch(filterlist,'*bok*') eq 0))
    filterlist = filterlist[keep]
    maggies = maggies[keep,*]
    ivarmaggies = ivarmaggies[keep,*]
    
    Bw = (where(strmatch(filterlist,'*Bw*',/fold)))[0]
    Bwtoss = where(redshift gt 1.0)

    ch3 = (where(strmatch(filterlist,'*ch3*')))[0]
    ch4 = (where(strmatch(filterlist,'*ch4*')))[0]
    ivarmaggies[Bw,Bwtoss] = 0.0
    ivarmaggies[ch3,*] = 0.0
    ivarmaggies[ch4,*] = 0.0
;   index = [8,9]
;   index = (where(sample.kmag gt 0.0))[1:2]

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       isedfit, paramfile, maggies, ivarmaggies, redshift, result, $
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
