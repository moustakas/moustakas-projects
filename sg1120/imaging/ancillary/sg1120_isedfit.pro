pro sg1120_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber
; jm09jul06nyu - written

    iopath = sg1120_path(/isedfit)
    paramfile = iopath+'sg1120_isedfit.par'

; --------------------------------------------------
; build the models
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=datapath

stop
    
; --------------------------------------------------
; read the merged GOODS catalog
    if n_elements(suffix) eq 0 then suffix = ''

    catfile = getenv('GOODS_LF_DIR')+'/data/lf_goods_catalog'+suffix+'.fits.gz'
    splog, 'Reading '+catfile
    goods = mrdfits(catfile,1)
    galaxy = string(goods.id,format='(I4.4)')
    
;   index = lindgen(n_elements(goods))
    index = where(goods.zfinal gt 2.9 and (goods.zfinal lt 3.0))
    lf_goods_to_maggies, goods, maggies, ivarmaggies

; --------------------------------------------------
; do the fitting
    if keyword_set(isedfit) then begin
; everything
       outprefix = 'BVIzJHKirac'
       isedfit, paramfile, maggies, ivarmaggies, goods.zfinal, result, $
         iopath=datapath, outprefix=outprefix, nminphot=nminphot, $
         clobber=clobber, debug=0, index=index
;; drop irac
;       outprefix = 'BVIzJHK'
;       ivarmaggies[7:9,*] = 0.0
;       isedfit, paramfile, maggies, ivarmaggies, goods.zfinal, result, $
;         iopath=datapath, outprefix=outprefix, nminphot=nminphot, $
;         clobber=clobber, debug=0, index=index
    endif

; --------------------------------------------------
; make QAplots
    if keyword_set(qaplot) then begin
; everything
       outprefix = 'BVIzJHKirac'
       isedfit_qaplot, paramfile, isedfit, iopath=iopath, galaxy=galaxy, $
         index=index, clobber=clobber, outprefix=outprefix
;; drop irac
;       outprefix = 'BVIzJHK'
;       isedfit_qaplot, paramfile, isedfit, iopath=iopath, $
;         galaxy=galaxy, index=index, clobber=clobber, outprefix=outprefix
    endif

stop    
    
return
end
