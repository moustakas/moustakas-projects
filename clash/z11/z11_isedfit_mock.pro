pro z11_isedfit_mock, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, noirac=noirac
; jm12nov13siena - fit mocked up z=11 candidate galaxy photometry that
; includes deeper Spitzer data

    datapath = clash_path(/z11)
    isedpath = datapath+'isedfit/'
    isedfit_sfhgrid_dir = datapath+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_sfhgrid.par'
    supergrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_supergrid.par'

    if n_elements(supergrid) eq 0 then supergrid = 1

; gather the photometry
    zlog = 0
    igm = 1
    prefix = 'z11_mock'
    zminmax = [10.6,10.8]
    nzz = 3
    
; do some trickery!    
    cat = rsex(clash_path(/z11)+'M0647JD_mock.cat')
    cat = struct_addtags(cat,replicate({z: 0.0},2))
;   cat.z = randomu(seed,nzz)*(zminmax[1]-zminmax[0])+zminmax[0]
    cat.z = 10.7

    super = yanny_readone(supergrid_paramfile)
    nsuper = n_elements(super)
    struct_print, super

    filters = z11_filterlist()
    paramfile = isedpath+prefix+'_paramfile.par'
    write_isedfit_paramfile, filters, prefix=prefix, minz=zminmax[0], $
      maxz=zminmax[1], nzz=nzz, igm=igm, isedpath=isedpath, clobber=clobber

;   paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
;     format='(i2.2)')+'_isedfit.par'
;   z11_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
;     nzz=nzz, zlog=zlog, igm=igm, super=super[gg], filters=filters

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, super=super, isedpath=isedpath, $
            clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; do the fitting!
       if keyword_set(isedfit) then begin
          z11_to_maggies, cat, maggies, ivarmaggies

          isedfit, paramfile, maggies, ivarmaggies, cat.z, isedpath=isedpath, $
            super=super, clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, outprefix=outprefix
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
;         index = [76,404]
          yrange = [31,21.5]
          xrange = [2000,70000] & 
          xlog = 1
          isedfit_qaplot, paramfile, result, super=super, isedpath=isedpath, $ ; galaxy='JD1+JD2+JD3', $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, xrange=xrange, yrange=yrange, xlog=xlog
       endif
    endfor
    
return
end
