pro megaspitzer_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, noirac=noirac, killerplot=killerplot
; jm11nov08ucsd - fit the z=9.6 galaxy, megaspitzer

    isedpath = clash_path(/megaspitzer)+'isedfit/'
    isedfit_sfhgrid_dir = clash_path(/megaspitzer)+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/megaspitzer/megaspitzer_sfhgrid.par'

    if n_elements(supergrid) eq 0 then supergrid = 1

; gather the photometry
    prefix = 'megaspitzer'
    refcat = read_megaspitzer()
    zminmax = [3.8,9.6]
    if keyword_set(killerplot) then begin
;      cat = mrdfits(clash_path(/megaspitzer)+'killerplot_simulations.fits',1)
;      outprefix = 'killerplot_simulations'
       cat = mrdfits(clash_path(/megaspitzer)+'killerplot.fits',1)
       outprefix = 'killerplot'
    endif else begin
       cat = refcat
       outprefix = 'alldropouts'
    endelse    
    if keyword_set(noirac) then outprefix = outprefix+'_noirac'

    nzz = 50
    zlog = 0
    igm = 1
    
    super = get_megaspitzer_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       megaspitzer_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
         nzz=nzz, zlog=zlog, igm=igm, super=super[gg], filters=filters
       
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; do the fitting!
       if keyword_set(isedfit) then begin
          if keyword_set(killerplot) then begin
             maggies = cat.maggies
             ivarmaggies = cat.ivarmaggies
          endif else begin
             megaspitzer_to_maggies, cat, maggies, ivarmaggies
          endelse
          
          if keyword_set(noirac) then begin
             isirac = where(strmatch(filters,'*irac*',/fold))
             ivarmaggies[isirac,*] = 0
          endif
          isedfit, paramfile, maggies, ivarmaggies, cat.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, outprefix=outprefix
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
;         index = lindgen(50)
          yrange = [31,21.5]
          xrange = [2000,70000] & xlog = 1
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=refcat.galaxy, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, xrange=xrange, yrange=yrange, xlog=xlog
       endif
    endfor
    
return
end
