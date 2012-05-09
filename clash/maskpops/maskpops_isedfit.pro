pro maskpops_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
  nzz=nzz, zlog=zlog, igm=igm, super=super, filters=filters
    
    filters = clash_filterlist()

    sfhgridstring = strtrim(super.sfhgrid,2)
    redcurvestring = redcurve2string(super.redcurve)
    synthmodels = strtrim(super.synthmodels,2)
    imf = strtrim(super.imf,2)
    h100 = clash_h100(omega0=omega0,omegal=omegal)
    
    splog, 'Writing '+paramfile
    zrange = string(zminmax[0],format='(F4.2)')+','+string(zminmax[1],$
      format='(F4.2)')+','+strtrim(nzz,2)+','+strtrim(zlog,2)+' # [minz,maxz,dz,log?]'
    openw, lun, paramfile, /get_lun
    printf, lun, 'h100                 '+string(h100,format='(F4.2)')
    printf, lun, 'omega0               '+string(omega0,format='(F4.2)')
    printf, lun, 'omegal               '+string(omegal,format='(F4.2)')
    printf, lun, 'synthmodels          '+synthmodels
    printf, lun, 'imf                  '+imf
    printf, lun, 'sfhgrid              '+sfhgridstring
    printf, lun, 'redcurve             '+redcurvestring
    printf, lun, 'prefix               '+prefix
    printf, lun, 'redshift             '+zrange
    printf, lun, 'igm                  '+strtrim(igm,2)+' # [0=no, 1=yes]'
    printf, lun, 'maxold               0 # [0=no, 1=yes]'
    printf, lun, 'filterlist           '+strjoin(filters,',')
    free_lun, lun 

return
end

pro maskpops_isedfit, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, noirac=noirac, lowz=lowz
; jm12may08ucsd

    isedpath = maskpops_path(/isedfit)
    isedfit_sfhgrid_dir = maskpops_path(/montegrids)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/maskpops/maskpops_sfhgrid.par'

; gather the photometry
    cat = read_maskpops()

    prefix = 'maskpops'
    zminmax = [fix(min(cat.z*10)),ceil(max(cat.z*10))]/10.0
    nzz = 3
    zlog = 0
    igm = 1
    
    super = get_maskpops_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

; loop on each supergrid
    for gg = 0, nsuper-1 do begin
       splog, 'working on grid '+strtrim(super[gg].supergrid,2)

       paramfile = isedpath+prefix+'_supergrid'+string(super[gg].supergrid,$
         format='(i2.2)')+'_isedfit.par'
       maskpops_write_paramfile, paramfile, prefix=prefix, zminmax=zminmax, $
         nzz=nzz, zlog=zlog, igm=igm, super=super[gg], filters=filters
       
; build the models
       if keyword_set(models) then begin
          isedfit_models, paramfile, iopath=isedpath, clobber=clobber, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif

; do the fitting!
       if keyword_set(isedfit) then begin
          maskpops_to_maggies, cat, maggies, ivarmaggies
          isedfit, paramfile, maggies, ivarmaggies, cat.z, iopath=isedpath, $
            clobber=clobber, sfhgrid_paramfile=sfhgrid_paramfile, $
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir
       endif       

; make some QAplots
       if keyword_set(qaplot) then begin
          yrange = [30,19]
;         xrange = [2000,70000] & xlog = 1
          isedfit_qaplot, paramfile, result, iopath=isedpath, galaxy=cat.prefix, $
            index=index, clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, $
            outprefix=outprefix, xrange=xrange, yrange=yrange, /xlog
       endif
    endfor
    
return
end
