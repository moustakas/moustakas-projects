pro macs0329_z6arcs_isedfit, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, noirac=noirac
; jm11nov08ucsd - fit the z~6 arcs in macs0329

    isedpath = clash_path(/ised)
    datapath = clash_path(/macs0329_z6arcs)

    isedfit_sfhgrid_dir = clash_path(/monte)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/clash_sfhgrid.par'

; read the supergrid parameter file
    supergrid = 3 ; current supergrids for this project
;   supergrid = [2,3] ; current supergrids for this project
    super = get_clash_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

; gather the photometry
    cat = read_z6arcs_photometry(adi=adi)

    igm = '1'
    nzz = '3'
    zlog = '0'
    prefix = 'macs0329_z6arcs'
    zminmax = [6.0,6.2]

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
          if keyword_set(noirac) then begin
             irac = where(strmatch(filters,'*irac*',/fold))
             cat.ivarmaggies[irac] = 0.0
          endif
          isedfit, paramfile, cat.maggies, cat.ivarmaggies, adi.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir;, index=index
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
          yrange = [32,21]
;         xrange = [2000,17000] & xlog = 0
          xrange = [2000,70000] & xlog = 1
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=adi.id, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, xrange=xrange, yrange=yrange, xlog=xlog
       endif
    endfor

return
end
