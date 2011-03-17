pro mz_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber, sdss=sdss, debug=debug
; jm10jul29ucsd - derive stellar masses for the MZ AGES and SDSS
; samples 

    if keyword_set(sdss) then prefix = 'sdss' else prefix = 'ages'
    iopath = ages_path(/projects)+'mz/isedfit/'
    paramfile = iopath+'mz_'+prefix+'_isedfit.par'
    params = read_isedfit_paramfile(paramfile)
    
; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; do the fitting!  
    if keyword_set(isedfit) then begin
       if keyword_set(sdss) then begin
          vagc = mz_get_vagc(sample=sample,letter=letter,poststr=poststr)
          post = read_vagc_garching(sample=sample,$
            letter=letter,poststr=poststr,/postlss)
          maggies = mz_get_maggies(post,/sdss,ivarmaggies=ivarmaggies)
          zobj = post.z
       endif else begin
          phot = read_ages(/photo)
          index = where((phot.imain eq 1) and (phot.z gt 0.01) and (phot.z lt 1.0))
;         splog, 'testing!!!!!'
;         index = index[900:949]
          maggies = mz_get_maggies(phot,ivarmaggies=ivarmaggies)
          zobj = phot.z
       endelse

       if (params.sfhgrid) lt 10 then galchunk = 300L else $
         galchunk = 100L

       isedfit, paramfile, maggies, ivarmaggies, zobj, result, $
         iopath=iopath, outprefix=outprefix, galchunksize=galchunk, $
         clobber=clobber, debug=debug, index=index
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       index = 500+lindgen(30)
;      splog, 'testing!!!!!'
;      phot = read_ages(/photo)
;      index = where((phot.imain eq 1) and (phot.z gt 0.01) and (phot.z lt 1.0))
;      index = index[900:949]
       isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
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
