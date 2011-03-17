pro sfrm_sdss_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber, debug=debug
; jm10feb16nyu - measure stellar masses for the SDSS galaxies

    common sfrm_sdss, lowz, galex, twomass
    
    iopath = ages_path(/project)+'sfrm/sdss_isedfit/'
    paramfile = iopath+'sdss_isedfit.par'

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; read the SDSS photometry - need to add 2MASS photometry!
    if (n_elements(lowz) eq 0) then begin
       lowzpath = sdss_path(/lowz)
       lowz = mrdfits(lowzpath+'lowz_catalog.dr6.fits.gz',1)
       galex = mrdfits(lowzpath+'lowz_galex_gr5.fits.gz',1)
;      galex = mrdfits(lowzpath+'lowz_galex.dr6.fits.gz',1)
       twomass = mrdfits(lowzpath+'lowz_twomass.dr6.fits.gz',1)
    endif
    ngal = n_elements(lowz)

    im_galex_to_maggies, galex, gmaggies, givarmaggies
;   galex_to_maggies, galex, gmaggies, givarmaggies
    sdss_to_maggies, smaggies, sivarmaggies, calib=lowz
    twomass_to_maggies, twomass, tmaggies, tivarmaggies
    
    maggies = [gmaggies,smaggies,tmaggies]
    ivarmaggies = [givarmaggies,sivarmaggies,tivarmaggies]
    redshift = lowz.zdist

;   index = where((maggies[1,*] gt 0.0))
;   index = where(total(maggies[7:9,*] gt 0.0,1) ge 1)
;   index = where((maggies[1,*] gt 0.0) and $
;     (total(maggies[7:9,*] gt 0.0,1) ge 1))

;; test
;    bmass = reform(k_sdss_bell(lowz[index].absmag[0:4])) ; stellar mass
;    bmass = alog10(bmass/0.7^2) ; h=1-->h=0.7
;    vmax = lowz[index].vmax/0.7^3 ; h=1-->h=0.7
;    binsize = 0.1
;    mf_vmax, bmass, 1.0/vmax/binsize, phi=phi, $
;      binmass=binmass, minmass=8.0
;    plot, binmass, phi, psym=10, xr=[7.5,12], yr=[1E-5,0.1], /ylog
;    oplot_bgd08_mf, maxis1, params=params, log=0, color='red', line=5
    
; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       isedfit, paramfile, maggies, ivarmaggies, redshift, result, $
         iopath=iopath, outprefix=outprefix, nminphot=nminphot, $
         clobber=clobber, debug=debug, index=index
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
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


;;    if keyword_set(doplots) then begin
;;       
;;       params = read_isedfit_paramfile(paramfile,iopath=iopath)
;;       fp = isedfit_filepaths(params,outprefix=outprefix,iopath=datapath)
;;       massrange = [6.5,13.0]
;;       splog, 'Reading '+fp.iopath+fp.isedfit_outfile
;;       res = mrdfits(fp.iopath+fp.isedfit_outfile+'.gz',1,/silent)
;;       
;;; K-correct
;;       djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=massrange, $
;;         yrange=massrange, xtitle='log (M/M_{\odot}) [K-correct]', $
;;         ytitle='log (M/M_{\odot}) [isedfit]'
;;       oploterror, ages[index].mass, res[index].mass50-0.25, res[index].mass_err, psym=8
;;;      djs_oplot, ages[index].mass, res.mass-0.25, psym=8, color='red'
;;;      djs_oplot, ages[index].mass, res.mass-0.25, psym=8
;;       djs_oplot, !x.crange, !y.crange, line=0, thick=2, color='red'
;;       im_legend, 'SFHGRID-'+string(params.sfhgrid,format='(I2.2)'), $
;;         /left, /top, box=0
;;       ss = im_stats(ages[index].mass-(res[index].mass50-0.25),/verbose)
;;       cc = get_kbrd(1)
;;       
;;       plot, res[index].age50, ages[index].mass-(res[index].mass50-0.25), ps=8, ysty=3, xsty=3
;;stop       
;;       
;;    endif
;;    
