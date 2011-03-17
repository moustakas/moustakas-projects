pro test_sfrm_sdss_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber, debug=debug
; jm10feb16nyu - measure stellar masses for the SDSS galaxies

    common sfrm_sdss, lowz, galex, twomass
    
    iopath = ages_path(/project)+'sfrm/sdss_isedfit/test/'
    paramfile = iopath+'sdss_isedfit.par'

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; read the SDSS photometry
    if (n_elements(lowz) eq 0) then begin
       lowz = mrdfits(sdss_path(/lowz)+'lowz_catalog.dr6.fits.gz',1)
       galex = mrdfits(sdss_path(/lowz)+'lowz_galex.dr6.fits.gz',1)
       twomass = mrdfits(sdss_path(/lowz)+'lowz_twomass.dr6.fits.gz',1)
    endif
    ngal = n_elements(lowz)

    bmass = reform(k_sdss_bell(lowz.absmag[0:4])) ; stellar mass
    bmass = alog10(bmass/0.7^2) ; h=1-->h=0.7
    vmax = lowz.vmax/0.7^3 ; h=1-->h=0.7
    binsize = 0.1
    mf_vmax, bmass, 1.0/vmax/binsize, phi=phi, $
      binmass=binmass, minmass=8.0
    plot, binmass, phi, psym=10, xr=[7.5,12], yr=[1E-5,0.1], /ylog
    oplot_bgd08_mf, maxis1, params=params, log=0, color='red', line=5
    index = where(bmass gt 11.4)
    
;    ised = mrdfits('ised.fits.gz',1)
;    rest = mrdfits('rest.fits.gz',1)
;    good = where((bmass gt 11.0) and (rest.maggies[1] gt 0) and $
;      (total(rest.maggies[7:9] gt 0.0,1) ge 1))
;
;    jv2ab = k_vega2ab(filterlist='twomass_J.par',/kurucz,/silent)
;    nuvmr = rest[good].galex_absmag[1]-rest[good].ugriz_absmag[2]
;    rmj = rest[good].ugriz_absmag[2]-(rest[good].ubvrijhk_absmag[5]+jv2ab)
;    qq = select_quiescent(nuvmr,rmj,active=aa)
;    index = good[qq]
;;   index = index[[6,15]]
;    plot, rmj, nuvmr, ps=6, sym=0.2, xr=[0,1.7], yr=[0.2,6.5]
    
;   index = where(bmass gt 10.9 and (bmass-ised.mass50 gt 0.4) and (ised.chi2 lt 50))
;   new = mrdfits('ugrizJHKs_chab_calzetti_sfhgrid02.fits.gz',1,rows=index)
    
    im_galex_to_maggies, galex, gmaggies, givarmaggies
    sdss_to_maggies, smaggies, sivarmaggies, calib=lowz
    twomass_to_maggies, twomass, tmaggies, tivarmaggies
    
    maggies = [gmaggies,smaggies,tmaggies]
    ivarmaggies = [givarmaggies,sivarmaggies,tivarmaggies]
    redshift = lowz.zdist

;   index = lindgen(ngal)
;   these = (where((maggies[0,index] gt 0.0) and $
;     (maggies[1,index] gt 0.0)))[3000:3020]
;   index = index[these]
    
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

    new = mrdfits('ugrizJHKs_bc03_chab_calzetti_sfhgrid02.fits.gz',1,rows=index)

    stop
    
    mf_vmax, new.mass_avg, 1.0/vmax[index]/binsize, phi=phi, binmass=binmass, minmass=8.0
    plot, binmass, phi, psym=10, xr=[7.5,12], yr=[1E-5,0.1], /ylog
    oplot_bgd08_mf, maxis1, params=params, log=0, color='red', line=5

stop    
    
    xr = [9.5,12]
    plot, new.mass50, bmass[index], ps=6, xr=xr, yr=xr & oplot, !x.crange, !y.crange
    plot, new.mass_avg, bmass[index], ps=6, xr=xr, yr=xr & oplot, !x.crange, !y.crange
    plot, new.ebv, new.mass-bmass[index], ps=6, xr=[0,0.4], yr=[-1,1]

stop    
    
    plot, new.mass50, bmass[index], ps=6, xr=xr, yr=xr & oplot, !x.crange, !y.crange
    plot, new.mass50, new.mass, ps=6, xr=xr, yr=xr & oplot, !x.crange, !y.crange
    plot, new.ebv_avg, new.mass_avg-bmass[index], ps=6, xr=[0,0.2], yr=[-2,2]

    
stop
    
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
