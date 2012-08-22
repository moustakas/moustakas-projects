pro z11_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, noirac=noirac, loz=loz
; jm12aug14siena - fit the z=11 candidate

    isedpath = clash_path(/z11)+'isedfit/'
    isedfit_sfhgrid_dir = clash_path(/z11)+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_sfhgrid.par'
    supergrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_supergrid.par'

    if n_elements(supergrid) eq 0 then supergrid = 1

; gather the photometry
    zlog = 0
    igm = 1
    if keyword_set(loz) then begin
       prefix = 'z11_loz'
       zminmax = [2.0,3.0]
       nzz = 10
    endif else begin
       prefix = 'z11_hiz'       
       zminmax = [10.0,12.0]
       nzz = 20
    endelse
    
; do some trickery!    
    cat = rsex(clash_path(/z11)+'M0647JD_final2sum.cat')
    nzz = 100
    cat = replicate(struct_addtags(cat[4],{z: 0.0}),nzz)
    cat.z = randomu(seed,nzz)*(zminmax[1]-zminmax[0])+zminmax[0]

    super = yanny_readone(supergrid_paramfile)
    nsuper = n_elements(super)
    struct_print, super

; while the computer is churning i have to (apologetically because i
; don't typically make this an issue, but i've gotten burned in the
; past) bring up the author order again; i only ask that i be placed
; in proportion to the amount of work i've done and/or its
; significance/importance for the paper; if that's N=15 so be it and
; i'm perfectly happy.  the reason i mention this is because i did
; *all* the metallicity work on what is now a highly-cited
; "metallicities of starburst galaxies" paper and i didn't speak up
; and stayed Nth author, so it's important that it
    
; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       z11_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
         nzz=nzz, zlog=zlog, igm=igm, super=super[gg], filters=filters
       
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; do the fitting!
       if keyword_set(isedfit) then begin
          z11_to_maggies, cat, maggies, ivarmaggies
          
          isedfit, paramfile, maggies, ivarmaggies, cat.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, outprefix=outprefix
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
;         index = [76,404]
          yrange = [31,21.5]
          xrange = [2000,70000] & 
          xlog = 1
          isedfit_qaplot, paramfile, result, iopath=isedpath, $ ; galaxy='JD1+JD2+JD3', $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, xrange=xrange, yrange=yrange, xlog=xlog
       endif
    endfor
    
return
end
