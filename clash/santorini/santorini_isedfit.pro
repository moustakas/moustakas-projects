pro santorini_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, noirac=noirac, lowz=lowz
; jm11nov08ucsd - fit the z=9.6 galaxy, santorini

    isedpath = clash_path(/santorini)+'isedfit/'
    isedfit_sfhgrid_dir = clash_path(/santorini)+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/santorini/santorini_sfhgrid.par'

; gather the photometry
    cat = read_santorini()

; consider both the z=9.6 and z=3.23 solutions
    if n_elements(supergrid) eq 0 then supergrid = [1,2,3,4,6]
    check = where(supergrid eq 5)
    if check[0] ne -1 and (keyword_set(lowz) eq 0) then $
      message, 'SUPERGRID=5 requires /LOWZ'
    if keyword_set(lowz) then begin
       prefix = 'santorini_lowz'
       cat.z = 3.2
       zminmax = [3.1,3.3]
       supergrid = 5
    endif else begin
       prefix = 'santorini'
       zminmax = [9.5,9.7]
       if n_elements(supergrid) eq 0 then supergrid = [1,2,3,4,6]
    endelse
    nzz = 3
    zlog = 0
    igm = 1
    
    super = get_santorini_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       clash_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
         nzz=nzz, zlog=zlog, igm=igm, super=super[gg], /useirac, $
         filters=filters
       
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; do the fitting!
       if keyword_set(isedfit) then begin
          santorini_to_maggies, cat, maggies, ivarmaggies
;         ivarmaggies[15] = 0
          isedfit, paramfile, maggies, ivarmaggies, cat.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
          yrange = [31,22]
          xrange = [2000,70000] & xlog = 1
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=cat.galaxy, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, xrange=xrange, yrange=yrange, xlog=xlog
       endif
    endfor
    
return
end
