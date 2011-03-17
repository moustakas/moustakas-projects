pro irclusters_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  chi2grid_qaplot=chi2grid_qaplot, measure=measure, clobber=clobber, $
  debug=debug
; jm10oct26ucsd - measure stellar masses for galaxies in
; Brodwin's sample of IR-selected clusters

    iopath = ages_path(/projects)+'irclusters/'
    paramfile = iopath+'irclusters_isedfit.par'
    param = read_isedfit_paramfile(paramfile)

    cat = read_irclusters(maggies=maggies,ivarmaggies=ivarmaggies)
;   test = lindgen(200)+500
;   cat = cat[test]
;   maggies = maggies[*,test]
;   ivarmaggies = ivarmaggies[*,test]
    ngal = n_elements(cat)

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; do the fitting!  do not use the NEWFIRM photometry nor IRAC/ch3-4
    if keyword_set(isedfit) then begin
       toss = where(strmatch(param.filterlist,'*ch3*') or $
         strmatch(param.filterlist,'*ch4*'))
       ivarmaggies[toss,*] = 0.0

;; fit at both the photometric redshift...
;       index = where((cat.z ge min(param.redshift)) and $
;         (cat.z le max(param.redshift)))
;       isedfit, paramfile, maggies, ivarmaggies, cat.z, $
;         iopath=iopath, clobber=clobber, index=index
;; ...and at the cluster redshift...
;;      res = mrdfits('BwRIJHKsirac_zclust_bc03_chab_calzetti_sfhgrid02.fits.gz',1)
;;      index = (where((res.mass_err lt 0.01) and (res.zobj gt 1.2) and (res.zobj lt 1.3) ))[0:49]
;       isedfit, paramfile, maggies, ivarmaggies, cat.zclust, $
;         iopath=iopath, outprefix=param.prefix+'_zclust', $
;         clobber=clobber;, index=index
;; ...and then refit at the photometric redshift excluding IRAC 
       index = where((cat.z ge min(param.redshift)) and $
         (cat.z le max(param.redshift)))
;      index = [200,239]
       alltoss = where(strmatch(param.filterlist,'*spitzer*'))
       ivarmaggies[alltoss,*] = 0.0
       isedfit, paramfile, maggies, ivarmaggies, cat.z, $
         iopath=iopath, outprefix=repstr(param.prefix,'irac',''), $
         clobber=clobber, index=index
    endif 

; --------------------------------------------------
; make some QAplots of a random subset of objects
    if keyword_set(qaplot) then begin
       qaindx = long(randomu(seed,100)*ngal)
       qaindx = qaindx[uniq(qaindx,sort(qaindx))]
       galaxy = 'Galaxy '+string(qaindx,format='(I4.4)')
       
; includes irac, photoz
       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
         galaxy=galaxy, clobber=clobber, outprefix=outprefix
; includes irac, zclust
       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
         galaxy=galaxy, clobber=clobber, outprefix=param.prefix+'_zclust'
; no irac, photoz
       isedfit_qaplot, paramfile, iopath=iopath, index=qaindx, $
         galaxy=galaxy, clobber=clobber, outprefix=repstr(param.prefix,'irac','')
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
