pro test_sfrm_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber, debug=debug
; jm10feb05nyu - measure stellar masses for the AGES/SFRM galaxies  

    iopath = ages_path(/projects)+'sfrm/isedfit/test/'
    paramfile = iopath+'sfrm_isedfit.par'

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath

; --------------------------------------------------
; read the AGES photometry
    maggies = read_ages_photometry(ivarmaggies=ivarmaggies,$
      index=index,filterlist=filterlist,redshift=redshift,$
      /usegalex)
;   these = (where((maggies[0,index] gt 0.0) and $
;     (maggies[1,index] gt 0.0)))[3000:3200]
;   index = index[these]

    ised = mrdfits('ised.fits.gz',1)
    rest = mrdfits('rest.fits.gz',1)
    index = where(ised.ebv50 gt 0.4 and (rest.ugriz_absmag[0]-rest.ugriz_absmag[1]) gt 1.8) 
    
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
;      these = (where((maggies[0,index] gt 0.0) and (maggies[1,index] gt 0.0)))[3000:3100]
;      index = index[these]
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
